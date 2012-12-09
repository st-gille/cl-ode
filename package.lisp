(defpackage :ode
  (:use :cl :cffi :alexandria)
  (:export :runge-kutta :newton-test))

(in-package :ode)

(defun nested (lst)
  (if (null (rest lst))
    (car lst)
    (list (car lst) (nested (rest lst)))))

(defun apply-predicate (object pred)
  (if (atom pred)
    (list pred object)
    (nested (append pred (list object)))))

(defmacro apply-helper (key object &rest predicates)
  (once-only (object)
    (let1 (tests (mapcar (curry #'apply-predicate object) predicates))
      `(,key ,@tests))))

(defmacro all-p (object &rest predicates)
  `(apply-helper and ,object ,@predicates))

(defmacro any-p (object &rest predicates)
  `(apply-helper or ,object ,@predicates))

(defmacro aif (test-form then-form &optional else-form)
  `(let ((it ,test-form))
        (if it ,then-form ,else-form)))

(defmacro let1 (binding &body body)
  `(let (,binding)
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

(defmacro make-alist-calling (fname &rest rest)
  "Construct an alist of symbols and arguments, calling the function fname on each argument."
  `(list ,@(mapcar (lambda (x) `(cons (quote ,(car x)) (,fname ,@(second x)))) rest)))

(defun modify-nth (lst new-val n)
  (cond
    ((null lst) (error "List too short."))
    ((< n 0) (error "Cannot access negative positions."))
    ((= n 0) (cons new-val (cdr lst)))
    (t (cons (car lst) (modify-nth (cdr lst) new-val (1- n))))))

(defun make-simple-list (dim &optional (formula #'+))
  (loop for j below dim collect (coerce (funcall formula j) 'double-float)))

(defun make-simple-matrix (rows cols &optional (formula #'+))
  (loop for i below rows collect (loop for j below cols collect (coerce (funcall formula i j) 'double-float))))

(defun list-to-2d-array (list)
  (make-array (list (length list)
                    (length (first list)))
              :initial-contents list))

(defun split (list n)
  (let ((res)
        (tmp))
    (loop for item in list
          for i from 1
          do (push item tmp)
          when (= n i)
          do (progn
                    (push (nreverse tmp) res)
                    (setq i 0)
                    (setf tmp nil)))
    (nreverse res)))
