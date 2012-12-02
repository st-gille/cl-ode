(in-package :ode)

(define-foreign-library libnewton
                        (:unix "libnewton.so.1.0.0")
                        (t (:default "libnewton")))
(use-foreign-library libnewton)

(defcfun "lr_decomp" :int8 (flags :uint64 (:pointer (:pointer :double)) (:pointer :uint64)))
(defcfun "lr_solve" :int8 (flags :uint64 (:pointer (:pointer :double)) (:pointer :double) (:pointer :double) (:pointer :uint64)))

(defmacro with-foreign-pointers (bindings &body body)
  "Wrap every binding (var size &optional size-var) in bindings with with-foreign-pointer."
  (if (null bindings)
    `(progn ,@body)
    (let ((inner-bind `(with-foreign-pointer ,(car bindings) ,@body)))
      (reduce (lambda (form binding) `(with-foreign-pointer ,binding ,form)) (rest bindings) :initial-value inner-bind))))

(with-foreign-pointer (y 255 strlen)
                      (with-foreign-pointer (x 5 dim)
                                            (setf (mem-ref x :double (1- dim)) 1.0d0)
                                            (setf (mem-ref y :double (1- strlen)) 1.0d0)))

(with-foreign-pointers ((x 4 dim) (y 255 strlen))
                       (setf (mem-ref x :double (1- dim)) 1.0d0)
                       (setf (mem-ref y :double (1- strlen)) 2.0d0))
