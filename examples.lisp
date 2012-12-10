(in-package :ode)

(defun identity-ode (tt x)
  x)

(defmacro test-runge-kutta (tableau &optional (method 'runge-kutta))
  `(with-tableau ,tableau
     (format t "Evaluating 100 steps with tableau ~A and ~A" ,tableau (quote ,method))
     (time (,method  #'identity-ode '(1.0d0 0.5d0 0.3d0) 0.00d0 1.0d0 0.01d0))
    t))

(defun sine-ode (tt xx)
  (list (second xx) (- (first xx))))

(defun caching-example ()
  (multiple-value-bind (sol graph)
    (cached-ode-solution #'sine-ode '(0.0 1.0) 0.0 0.1)
    (flet ((print-len ()
             (format t "Length of known graph: ~d"
                     (length (funcall graph))))
           (print-eval (at)
             (progn
               (print "Evaluating...")
               (format t "Value of solution at ~,5e: ~,4e~%"
                       at
                       (time (funcall sol at))))))
      (format t "Initital graph: ~A."
              (funcall graph))
      (with-tableau :explicit-heun
        (print-eval pi)
        ;(format t "Known graph: ~A~%" (funcall graph))
        (print-len)
        (print-eval (/ pi 2))
        (print-eval 1.0)
        ;(format t "Known graph: ~A~%" (funcall graph))
        (print-len)
        (print-eval (* 2 pi))
        (print-len))
      (print "Returning solution as primary, graph as secondary value.")
      (values sol graph))))

(defun run-examples ()
  (progn
    (test-runge-kutta :explicit-euler)
    (test-runge-kutta :explicit-euler runge-kutta-values)
    (test-runge-kutta :explicit-euler runge-kutta-graph)
    (test-runge-kutta :implicit-euler)
    (test-runge-kutta :implicit-trapezoid runge-kutta)
    (test-runge-kutta :classic)
    (caching-example)))
