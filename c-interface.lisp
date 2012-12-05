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

(defun set-c-vector (vec vals &optional (type :double))
  "Convert list or vector to c-array"
  (loop for val being the elements of vals
        for i from 0
        do (setf (mem-aref vec type i) val)))

(defun convert-from-c-vector (size vec &optional (type :double))
  (loop for i below size
        collect (mem-aref vec type i)))

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

(defmacro with-matrix ((name &key (type :double) initial-contents (dims nil have-dims)) &body body)
  (with-gensyms (g-dims g-contents g-rows g-cols)
    `(let* (,@(if have-dims `((,g-dims ,dims)))
            ,@(if (not have-dims) `((,g-contents ,initial-contents)))
            (,g-rows ,(if have-dims `(first ,g-dims) `(length ,g-contents)))
            (,g-cols ,(if have-dims `(second ,g-dims) `(length (first ,g-contents)))))
    (with-foreign-object (,name :pointer ,g-rows)
     (alloc-c-matrix ,name ,g-rows ,g-cols ,@(if (not have-dims) `(:initial-contents ,g-contents)) :type ,type)
     ,@body
     (free-c-matrix ,name ,g-rows)))))

(defcallback cb-f :void ((py :pointer) (px :pointer) (pres :pointer))
  (let ((y (mem-ref py :double))
        (x (mem-ref px :double)))
    (setf (mem-ref pres :double) (- 1 (* x (exp y))))))

(defcallback df :void ((py :pointer) (px :pointer) (ppres ::pointer))
  (let ((y (mem-ref py :double))
        (x (mem-ref px :double))
        (pres (mem-ref ppres :pointer)))
    (setf (mem-ref pres :double) (- (* x (exp y))))))

(defun newton-test ()
  (with-foreign-objects ((y :double) (x :double))
    (set-c-vector y '(0.6d0))
    (set-c-vector x '(1.5d0))
    (newtons-method 1 y x (get-callback 'cb-f) (get-callback 'df))
    (print-vector 1 y)))

(defun lr-test (rows)
  (with-foreign-objects ((x :double rows) (b :double rows) (pivot :uint64 rows))
    (with-matrix (A :initial-contents (make-simple-matrix rows rows (lambda (i j)(declare (ignore i j)) (random 100))))
      (print-matrix rows rows A)
      (if (= 1 (lr-decomp rows A pivot)) (print "unheil!\n")
        (progn
          (print-matrix rows rows A)
          (print-vector-u rows pivot )
          (set-c-vector b (make-simple-list rows (lambda (i)(declare (ignore i)) (random 100))))
          (print-vector rows b)
          (print-vector rows x)
          (lr-solve rows A b x pivot)
          (print-vector rows x))))))
