(in-package :ode)

(defun identity-ode (tt x)
  x)

(defmacro test-runge-kutta (tableau &optional (method 'runge-kutta))
  `(with-tableau ,tableau
     (time (,method  #'identity-ode '(1.0d0 0.5d0 0.3d0) 0.00d0 1.0d0 0.01d0))
    ))


(defun run-examples ()
  (test-runge-kutta explicit-euler))
;(test-runge-kutta explicit-euler runge-kutta-values)
;(test-runge-kutta explicit-euler runge-kutta-graph)
;(test-runge-kutta implicit-euler)
;(test-runge-kutta implicit-trapezoid runge-kutta-graph)
;(test-runge-kutta classic)
