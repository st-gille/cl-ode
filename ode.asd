(defsystem "ode"
 :description "Solve explicit, first order ODEs using Runge-Kutta methods."
 :version "0.1"
 :author "Stefan Gille <gille@numasoft.de>"
 :licence "Public Domain"
 :components ((:file "package")
              (:file "runge-kutta")
              (:file "c-interface")))
