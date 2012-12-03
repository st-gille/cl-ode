(require :cffi)
(defpackage :ode
  (:use :cl :cffi)
  (:export :runge-kutta))
