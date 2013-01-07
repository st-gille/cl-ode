(defpackage :ode
  (:use :cl :cl-user :cffi :alexandria)
  (:export :butcher-tableau
           :show-tableaus
           :with-tableau
           :runge-kutta-graph
           :runge-kutta
           :runge-kutta-standalone
           :runge-kutta-values
           :runge-kutta-values-standalone
           :cached-ode-solution
           :with-matrix
           :newton-solver
           :run-examples
           :caching-example))
