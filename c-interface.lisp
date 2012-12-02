(in-package :ode)

(define-foreign-library libnewton
                        (:unix "libnewton.so.1.0.0")
                        (t (:default "libnewton")))
(use-foreign-library libnewton)

(defcfun "lr_decomp" :int8 (flags :uint64 (:pointer (:pointer :double)) (:pointer :uint64)))
(defcfun "lr_solve" :int8 (flags :uint64 (:pointer (:pointer :double)) (:pointer :double) (:pointer :double) (:pointer :uint64)))
