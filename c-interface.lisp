(in-package :ode)

(define-foreign-library libnewton
                        (:unix "libnewton.so.1.0.0")
                        (t (:default "libnewton")))
(use-foreign-library libnewton)

(defcfun "lr_decomp" :int8
         (dim :uint64)
         (matrix (:pointer (:pointer :double)))
         (pivot (:pointer :uint64)))

(defcfun "lr_solve" :int8
         (dim :uint64)
         (decomp-matrix (:pointer (:pointer :double)))
         (rhs (:pointer :double))
         (solution (:pointer :double))
         (pivot (:pointer :uint64)))

(defcfun "newtons_method" :int8
         (dim :uint64)
         (x (:pointer :double))
         (y (:pointer :double))
         (func :pointer)
         (diff_func :pointer))

(defcfun "print_vector" :void (dim :uint64) (x (:pointer :double)))
(defcfun "print_matrix" :void
         (rows :uint64)
         (cols :uint64)
         (x (:pointer (:pointer :double))))

(defmacro wrap-binds (key bindings &body body)
  "Wrap every binding in bindings with keyword."
  (if (null bindings)
    `(progn ,@body)
    (let ((inner-bind `(,key ,(car bindings) ,@body)))
      (reduce (lambda (form binding) `(,key ,binding ,form)) (rest bindings) :initial-value inner-bind))))

(defun wrap-binds-test ()
 (macroexpand-1 '(wrap-binds with-foreign-pointer ((x dim) (b dim)) form1 form2)))

(defun set-c-vector (vec vals &optional (type :double))
  "Convert list or vector to c-array"
  (loop for val being the elements of vals
        for i from 0
        do (setf (mem-aref vec type i) val)))

(defun make-simple-list (dim &optional (formula #'+))
  (loop for j below dim collect (coerce (funcall formula j) 'double-float)))

(defun set-c-vector-test ()
  (let ((dim 5))
    (with-foreign-object (x :double dim)
      (print-vector dim x)
      (set-c-vector x (make-simple-list dim))
      (print-vector dim x)
      (set-c-vector x (make-array dim :initial-contents (make-simple-list dim (lambda (i) (+ 2 (* i i))))))
      (print-vector dim x))))
(set-c-vector-test)

(defmethod set-c-matrix (A rows)
  "Row-wise initialization of c-style matrices."
  (loop for row being the elements of rows
        for i from 0
        do (set-c-vector (mem-aref A :pointer i) row)))

(defmethod set-c-matrix (A (rows array))
  "Row-wise initialization of c-style matrices via 2d-arrays."
      (loop for i below (array-dimension rows 0)
            do (loop for j below (array-dimension rows 1)
                   do (setf (mem-aref (mem-aref A :pointer i) :double j) (aref rows i j)))))

(defun alloc-c-matrix (A rows cols &key initial-contents (type :double))
  "Allocate enough memory to hold <rows>*<cols> objects of foreign type <type>. Row-wise initialization possible."
  (loop for i below (if rows rows (length initial-contents))
        do
        (prin1 "oho")
        (setf (mem-aref A :pointer i) (if initial-contents
                                        (foreign-alloc type :initial-contents (elt initial-contents i))
                                        (foreign-alloc type :count cols)))))

(defun free-c-matrix (A rows)
  (dotimes (row rows)
    (foreign-free (mem-aref A :pointer row))))

(defun make-simple-matrix (rows cols &optional (formula #'+))
  (loop for i below rows collect (loop for j below cols collect (coerce (funcall formula i j) 'double-float))))

(defun list-to-2d-array (list)
  (make-array (list (length list)
                    (length (first list)))
              :initial-contents list))

(defun set-c-matrix-test ()
  (let* ((rows 3)
         (cols 4)
         (m1 (make-simple-matrix rows cols (lambda (i j) (random 100))))
         (m2 (list-to-2d-array (make-simple-matrix rows cols (lambda (i j) (+ i (* i j)))))))
    (with-foreign-object (A :double rows)
      (alloc-c-matrix A rows cols :initial-contents (make-simple-matrix rows cols))
      (print-matrix rows cols A)
      (set-c-matrix A m1)
      (print-matrix rows cols A)
      (set-c-matrix A m2)
      (print-matrix rows cols A)
      (free-c-matrix A rows))))

(set-c-matrix-test)

(defun make-let-env (binds &rest rest)
  `(let* ,(loop for (var bind) in binds
                 collect `(,var ,bind))
     ,@rest))

(defmacro with-matrix ((name &key (type :double) initial-contents (dims '(0 0) have-dims)) &body body)
  (let* ((g-rows (if have-dims (first dims) (gensym)))
         (g-cols (if have-dims (second dims) (gensym)))
         (g-contents (if have-dims initial-contents (gensym)))
         (rowform (if have-dims
                    `(first ,dims)
                    `(length ,g-contents)))
         (colform (if have-dims
                    `(second ,dims)
                    `(length (first ,g-contents))))
         (full-body `(with-foreign-object (,name :pointer ,g-rows)
                       (alloc-c-matrix ,name ,g-rows ,g-cols :initial-contents ,g-contents :type ,type)
                       ,@body
                       (free-c-matrix ,name ,g-rows))))
    (if have-dims
      full-body
      (make-let-env `((,g-contents ,initial-contents)
                      (,g-rows ,rowform)
                      (,g-cols ,colform))
                    full-body))))
