(in-package :ode)

(defvar *eps* 1.0e-14)
(defvar *step* 1.0e-4)

(defun central-diff-quot (f)
  "Return new function approximating the derivative of input."
  (lambda (x)
    (let* ((x+ (+ x *step*))
           (x- (- x *step*))
           (fx+ (funcall f x+))
           (fx- (funcall f x-))
           (formula (lambda (fx fy) (* 0.5 (/ 1.0 *step*) (- fx fy)))))
      (mapcar formula fx+ fx-))))

(defun make-autonomous (f)
    (lambda (x) (list 1 (funcall f (first x) (rest x)))))

(defstruct (butcher-tableau
  (:constructor make-butcher (&optional (matrix '((0))) (weights '(1)) (time-coeffs '(0)) (is-implicit nil))))
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
  (let ((m (if (symbolp tableau) `(cdr (assoc (quote ,tableau) *tableaus*)) tableau)))
    `(let ((*selected-tableau* ,m)
           (*A* ) (*b*) (*c*) (*is-implicit*))
       (update-variables *selected-tableau*)
       ,@body)))

(defun eval-point (x k a)
  (if k
    (mapcar (lambda (x k a)(+ x (* k a))) x k a)
    x))

(defun solve-explicit-rks (f x0 t0 stepsize)
  (let* ((k nil)
         (scale (lambda (x) (* stepsize x)))
         (substitution-step  (lambda (a c)
                               (mapcar scale
                                       (funcall f
                                                (+ t0 (* c stepsize))
                                                (eval-point x0 k a))))))
    (nconc k (mapcar substitution-step *A* *c*))))

(defun runge-kutta-single-step (f x0 t0 stepsize)
  (let* ((ks (solve-explicit-rks f x0 t0 stepsize))
         (bks (mapcar (lambda (bi k) (mapcar (lambda (ki) (* bi ki)) k)) *b* ks)))
    (apply #'mapcar (cons #'+ (cons x0 bks)))))

(defun make-discretization (t0 t1 stepsize)
  (loop as i from t0 to t1 by stepsize
        collect i))

(defun runge-kutta-graph (f x0 t0 t1 stepsize)
  (loop as next-t from t0 to t1 by stepsize
        for x = x0
        then (runge-kutta-single-step f x next-t stepsize)
        collect (list next-t x) into graph
        finally (return
                  (nconc graph
                         (list (list t1
                                     (let ((last-t (- next-t stepsize)))
                                       (runge-kutta-single-step f x last-t (- t1 last-t)))))))))

(defun runge-kutta (f x0 t0 t1 stepsize)
  (if (>= t0 t1)
    (list x0)
    (let ((next (runge-kutta-single-step f x0 t0 stepsize)))
      (cons x0 (runge-kutta f next (+ t0 stepsize) t1 stepsize)))))

(defun eval-f (f t0 t1 stepsize)
  (let ((ft0 (funcall f t0)))
    (if (> t0 t1)
      (list ft0)
      (cons ft0 (eval-f f (+ t0 stepsize) t1 stepsize)))))

(defun num-equal (a b)
  (and a b (< (abs (- a b)) *eps*)))

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
