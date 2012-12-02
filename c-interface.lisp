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
