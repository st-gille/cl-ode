(asdf:defsystem ode
 :depends-on (cffi)
 :description "Solve explicit, first order ODEs using Runge-Kutta methods."
 :version "0.1"
 :author "Stefan Gille <gille@numasoft.de>"
 :licence "Public Domain"
 :serial t
 :components ((:file "package")
              (:file "common")
              (:file "c-interface")
              (:file "runge-kutta")
              (:file "examples")))

