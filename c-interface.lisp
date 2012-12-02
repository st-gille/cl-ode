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

