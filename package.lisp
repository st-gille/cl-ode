(defpackage :ode
  (:use :cl :cffi)
  (:export :runge-kutta :newton-test))

(in-package :ode)

(defmacro with-gensyms ((&rest names) &body body)
  `(let ,(loop for n in names collect `(,n (gensym)))
     ,@body))

(defmacro with-gensym-if-not (cond (&rest bindings) &body body)
  (if cond
    `(let ,(loop for (var init-form) in bindings collect `(,var (,init-form))) ,@body)
    `(with-gensyms ,(mapcar #'car bindings) ,@body)))

(defmacro wrap-binds (form bindings &body body)
  "Nest all bindings in <bindings> within <form>."
  (if bindings
    `(,form ,(first bindings)
            (wrap-binds ,form ,(rest bindings) ,@body))
    `(progn ,@body)))

(defun wrap-binds-test ()
 (macroexpand '(wrap-binds with-something ((var1 init1) (var2 init2)) form1 form2)))

(defun make-simple-list (dim &optional (formula #'+))
  (loop for j below dim collect (coerce (funcall formula j) 'double-float)))

(defun make-simple-matrix (rows cols &optional (formula #'+))
  (loop for i below rows collect (loop for j below cols collect (coerce (funcall formula i j) 'double-float))))

(defun list-to-2d-array (list)
  (make-array (list (length list)
                    (length (first list)))
              :initial-contents list))
