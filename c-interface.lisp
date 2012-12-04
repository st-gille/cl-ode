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
(defcfun "print_vector_i" :void (dim :uint64) (x (:pointer :int)))
(defcfun "print_vector_u" :void (dim :uint64) (x (:pointer :uint64)))
(defcfun "print_matrix" :void
         (rows :uint64)
         (cols :uint64)
         (x (:pointer (:pointer :double))))

(defmacro wrap-binds (form bindings &body body)
  "Nest all bindings in <bindings> within <form>."
  (if bindings
    `(,form ,(first bindings)
            (wrap-binds ,form ,(rest bindings) ,@body))
    `(progn ,@body)))

(defun wrap-binds-test ()
 (macroexpand '(wrap-binds with-foreign-pointer ((x dim) (b dim)) form1 form2)))

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
        (setf (mem-aref A :pointer i) (if initial-contents
                                        (foreign-alloc type :initial-contents (elt initial-contents i))
                                        (foreign-alloc type :count cols)))))

(defun free-c-matrix (A no-of-rows)
  (dotimes (row no-of-rows)
    (foreign-free (mem-aref A :pointer row))))

(defun make-simple-matrix (rows cols &optional (formula #'+))
  (loop for i below rows collect (loop for j below cols collect (coerce (funcall formula i j) 'double-float))))

(defun list-to-2d-array (list)
  (make-array (list (length list)
                    (length (first list)))
              :initial-contents list))


(defmacro with-matrix ((name &key (type :double) initial-contents (dims '(0 0) have-dims)) &body body)
  (let* ((g-rows (if have-dims (first dims) (gensym)))
         (g-cols (if have-dims (second dims) (gensym)))
         (g-contents (if have-dims initial-contents (gensym)))
         (rowform (if (not have-dims)
                    `(length ,g-contents)))
         (colform (if (not have-dims)
                    `(length (first ,g-contents))))
         (full-body `(with-foreign-object (,name :pointer ,g-rows)
                       (alloc-c-matrix ,name ,g-rows ,g-cols :initial-contents ,g-contents :type ,type)
                       ,@body
                       (free-c-matrix ,name ,g-rows))))
    (if have-dims
      full-body
      `(let* ((,g-contents ,initial-contents)
              (,g-rows ,rowform)
              (,g-cols ,colform))
         ,full-body))))

(defun lr-test (rows)
  (with-foreign-objects ((x :double rows) (b :double rows) (pivot :uint64 rows))
    (with-matrix (A :initial-contents (make-simple-matrix rows rows (lambda (i j)(declare (ignore i j)) (random 100))))
      (print-matrix rows rows A)
      (when (= 1 (lr-decomp rows A pivot)) (print "unheil!\n")
        (print-matrix rows rows A)
        (print-vector-u rows pivot )
        (set-c-vector b (make-simple-list rows (lambda (i)(declare (ignore i)) (random 100))))
        (print-vector rows b)
        (print-vector rows x)
        (lr-solve rows A b x pivot)
        (print-vector rows x)))))
