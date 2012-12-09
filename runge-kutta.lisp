(in-package :ode)

(defun num-equal (a b)
  (and a b (< (abs (- a b)) *eps*)))

(defvar *eps* 1.0e-14)
(defvar *step* 1.0e-4)

(defun make-central-diff-quot (f)
  "Return new function approximating the derivative of input."
  (lambda (x) (central-diff-quot f x)))

(defun central-diff-quot (f x)
  (let* ((x+ (+ x *step*))
         (x- (- x *step*))
         (fx+ (funcall f x+))
         (fx- (funcall f x-))
         (formula (lambda (fx fy)  (* 0.5 (/ 1.0 *step*) (- fx fy)))))
    (if (atom fx-)
      (funcall formula fx+ fx-)
      (mapcar formula fx+ fx-))))

(defun make-jacobian (f)
  (lambda (x) (jacobian f x)))

(defun swap-matrix-layout (matrix)
  (apply #'mapcar (cons #'list matrix)))

(defun jacobian (f x)
  (let* ((dim (length x))
         (df-by-cols (loop for i below dim
                           collect (flet ((curry-at-nth (y)
                                                        (let1 (tmp (modify-nth x y i))
                                                          (funcall f tmp))))
                                     (central-diff-quot #'curry-at-nth (nth i x))))))
    (swap-matrix-layout df-by-cols)))

(defun make-autonomous (f)
  (lambda (x) (list 1 (funcall f (first x) (rest x)))))

(defstruct (butcher-tableau
             (:constructor make-butcher (&optional (matrix '((0)))
                                                   (weights '(1))
                                                   (time-coeffs '(0))
                                                   (is-implicit nil))))
  (matrix :read-only)
  (weights :read-only)
  (time-coeffs :read-only)
  (is-implicit :read-only))

(defparameter *tableaus*
  (make-alist-calling make-butcher
    (explicit-euler)
    (implicit-euler ('((1)) '(1) '(1) t))
    (explicit-heun ('((0) (1 0)) '(0.5 0.5) '(0 1)))
    (classic ('(() (0.5) (0 0.5) (0 0 1)) '(1/6 1/3 1/3 1/6) '(0 0.5 0.5 1)))
    (implicit-trapezoid ('((0 0) (0.5 0.5)) '(0.5 0.5) '(0 1) t))))

(defvar *selected-tableau*)
(defvar *A*)
(defvar *b*)
(defvar *c*)
(defvar *is-implicit*)

(defun update-variables (selected-tableau)
  (declare (type butcher-tableau selected-tableau))
  (with-slots (matrix weights time-coeffs is-implicit) selected-tableau
    (setf *A* matrix
          *b* weights
          *c* time-coeffs
          *is-implicit* is-implicit)))

(defmacro with-runge-kutta (tableau &body body)
  (let1 (m (if (symbolp tableau) `(cdr (assoc (quote ,tableau) *tableaus*)) tableau))
    `(let ((*selected-tableau* ,m)
           (*A* ) (*b*) (*c*) (*is-implicit*))
       (update-variables *selected-tableau*)
       ,@body)))

(defun dot-product (u v)
  (apply #'+ (mapcar #'* u v)))

(defun eval-point (x k aj)
  (if k
    (flet ((fma (xi ki) (+ xi (dot-product aj ki))))
      (apply #'mapcar (list #'fma x (swap-matrix-layout k))))
    x))

(defun scaled-ode-eval (ode stepsize ai ci k t0 x0)
  (mapcar (lambda (fi) (* stepsize fi))
          (funcall ode
                   (+ t0 (* ci stepsize))
                   (eval-point x0 k ai))))

(defun solve-explicit-rks (ode x0 t0 stepsize)
  (let ((k))
    (loop for ai in *A*
          for ci in *c*
          do (let ((ki (list (scaled-ode-eval ode stepsize ai ci k t0 x0))))
               (appendf k ki)))
    k))

(defun solve-implicit-rks (ode x0 t0 stepsize)
  (let* ((dim (length x0))
         (rk-steps  (length *c*))
         (len (* dim rk-steps)))
    (labels ((weight-diff (kij fj) (- kij (* fj stepsize)))
             (eqs-for-ki (ki ode-system) (mapcar #'weight-diff ki ode-system))
             (ode-system-for-step (k) (flet ((ode-eval (ai ci)
                                                       (funcall ode
                                                                (+ t0 (* ci stepsize))
                                                                (eval-point x0 k ai))))
                                        (mapcar #'ode-eval *A* *c*)))
             (implicit-system (k-flat) (let ((k (split k-flat dim)))
                                         (mapcan #'eqs-for-ki k (ode-system-for-step k)))))
      (split (newton-solver #'implicit-system len) dim))))

(defun runge-kutta-single-step (f x0 t0 stepsize)
  (let* ((solver (if *is-implicit* #'solve-implicit-rks #'solve-explicit-rks))
         (k (funcall solver f x0 t0 stepsize))
         (bk (mapcar (lambda (bi ki) (mapcar (lambda (kij) (* bi kij)) ki)) *b* k)))
    (apply #'mapcar (cons #'+ (cons x0 bk)))))

(defun ensure-last-step (ode last-x last-t t1)
  (if (or (>= last-t t1) (num-equal last-t t1))
    nil
    (runge-kutta-single-step ode last-x last-t (- t1 last-t))))

(defun runge-kutta-graph (ode x0 t0 t1 stepsize)
  (loop as next-t from t0 to t1 by stepsize
        for x = x0
        then (runge-kutta-single-step ode x next-t stepsize)
        collect (list next-t x) into graph
        finally (return
                  (nconc graph (aif (ensure-last-step ode x (- next-t stepsize) t1)
                                 (list (list t1 it)))))))

(defun runge-kutta-values (ode x0 t0 t1 stepsize)
  (mapcar #'second (runge-kutta-graph ode x0 t0 t1 stepsize)))

(defun runge-kutta-values-standalone (ode x0 t0 t1 stepsize)
  (loop as next-t from t0 to t1 by stepsize
        for x = x0
        then (runge-kutta-single-step ode x next-t stepsize)
        collect x into vals
        finally (return (nconc vals (aif (ensure-last-step ode x (- next-t stepsize) t1)
                                      (list it))))))

(defun runge-kutta (ode x0 t0 t1 stepsize)
  (cadar (last (runge-kutta-graph ode x0 t0 t1 stepsize))))

(defun runge-kutta-standalone (ode x0 t0 t1 stepsize)
  (loop as next-t from t0 to t1 by stepsize
        for x = x0
        then (runge-kutta-single-step ode x next-t stepsize)
        finally (return (aif (ensure-last-step ode x (- next-t stepsize) t1)
                          it
                          x))))

(defun make-graph (f t0 t1 stepsize)
  (let ((ft0 (funcall f t0)))
    (if (> t0 t1)
      (list ft0)
      (cons ft0 (eval-f f (+ t0 stepsize) t1 stepsize)))))

(defun nassoc (value alist)
  (assoc value alist :test #'num-equal))

(defun find-last-smaller (item alist)
  (loop for i = alist
        then (rest i)
        when (or (null (rest i))
                 (> (caar (rest i)) item))
        return i))

(defun insert-at (place-to-put new-points)
  (let ((new-tail (cdr place-to-put)))
    (if new-tail
      (setf (cdr place-to-put) (append new-points new-tail))
      (nconc place-to-put new-points))))

(defun cached-ode-solution (the-ode init-t0 init-x0 init-stepsize)
  (let ((stepsize (max *eps* init-stepsize))
        (known-graph (list (list init-t0 init-x0))))
    (values
      (lambda (t1)
        (let* ((closest (find-last-smaller t1 known-graph))
               (t0 (caar closest))
               (x0 (cadar closest)))
          (if (< t1 t0)
            (error "Can't go back in time."))
          (cond
            ((num-equal t0 t1)
             x0)
            (t
              (let ((graph (runge-kutta-graph the-ode x0 t0 t1 stepsize)))
                (insert-at closest graph)
                (cadar (last graph)))))))
      (lambda () known-graph))))
