(in-package :ode)

(defstruct (butcher-tableau
             (:constructor make-butcher
                           (&optional (matrix '((0)))
                                      (weights '(1))
                                      (time-coeffs '(0))
                                      (is-implicit (not (lower-triangular-p matrix))))))
  (matrix :read-only)
  (weights :read-only)
  (time-coeffs :read-only)
  (is-implicit :read-only))

(defparameter *tableaus*
  (make-alist-calling make-butcher
    (:explicit-euler)
    (:implicit-euler
      ('((1)) '(1) '(1) t))
    (:explicit-heun
      ('(() (1 0)) '(0.5 0.5) '(0 1) nil))
    (:classic
      ('(() (0.5) (0 0.5) (0 0 1)) '(1/6 1/3 1/3 1/6) '(0 0.5 0.5 1) nil))
    (:implicit-trapezoid
      ('((0 0) (0.5 0.5)) '(0.5 0.5) '(0 1) t))))

;define known tableau names as keywords
(mapcar #'make-keyword (mapcar #'string (mapcar #'first *tableaus*)))

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

(update-variables (make-butcher))

(defmacro with-tableau (tableau &body body)
  "Create context that defines a specific butcher tableau.
  Pass either an object of type butcher-tableau or a symbol.
  If passed a symbol, it's associated tableau will
  be retrieved from a pre-defined list."
  (let1 (m (if (symbolp tableau)
             `(cdr (assoc (quote ,tableau) *tableaus*))
             tableau))
    `(let1 (*selected-tableau* ,m)
       (if (any-p *selected-tableau* null (not butcher-tableau-p))
         (error "Not a valid butcher-tableau: ~A~%" *selected-tableau*))
       (let ((*A* ) (*b*) (*c*) (*is-implicit*))
         (update-variables *selected-tableau*)
         ,@body))))

(defun eval-point (x k aj)
  (if k
    (flet ((fma (xi ki) (+ xi (dot-product aj ki))))
      (apply #'mapcar (list #'fma x (swap-matrix-layout k))))
    x))

(defun ode-eval (ode stepsize k x0 t0 ai ci)
  (funcall ode
           (+ t0 (* ci stepsize))
           (eval-point x0 k ai)))

(defun scaled-ode-eval (ode stepsize k x0 t0 ai ci)
  (mapcar (lambda (fi) (* stepsize fi))
          (ode-eval ode stepsize k x0 t0 ai ci)))

(defun solve-explicit-rks (ode x0 t0 stepsize)
  (let1 (k)
    (loop for ai in *A*
          for ci in *c*
          do (let1 (ki (list (scaled-ode-eval ode stepsize k x0 t0 ai ci)))
               (appendf k ki)))
    k))

(defun solve-implicit-rks (ode x0 t0 stepsize)
  (let* ((dim (length x0))
         (rk-steps  (length *c*))
         (len (* dim rk-steps)))
    (labels ((weight-diff (kij fj)
               (- kij (* fj stepsize)))
             (eqs-for-ki (ki ode-system)
               (mapcar #'weight-diff
                       ki
                       ode-system))
             (ode-system-for-step (k)
               (mapcar (curry #'ode-eval ode stepsize k x0 t0)
                       *A*
                       *c*))
             (implicit-system (k-flat)
               (let1 (k (split-n dim k-flat))
                 (mapcan #'eqs-for-ki
                         k
                         (ode-system-for-step k)))))
      (split-n dim (newton-solver #'implicit-system len)))))

(defun runge-kutta-single-step (ode x0 t0 stepsize)
  "Solve a single step of the ode with the current butcher tableau."
  (let* ((solver (if *is-implicit* #'solve-implicit-rks #'solve-explicit-rks))
         (k (funcall solver ode x0 t0 stepsize))
         (bk (mapcar (lambda (bi ki)
                       (mapcar (lambda (kij) (* bi kij))
                               ki))
                     *b*
                     k)))
    (apply #'mapcar (cons #'+ (cons x0 bk)))))

(defun ensure-last-step (ode last-x last-t t1)
  (if (or (>= last-t t1) (num-equal last-t t1))
    nil
    (runge-kutta-single-step ode last-x last-t (- t1 last-t))))

(defun runge-kutta-graph (ode x0 t0 t1 stepsize)
  "Calculate the graph of the solution of the initial value problem."
  (loop as next-t from t0 to t1 by stepsize
        for x = x0
        then (runge-kutta-single-step ode x next-t stepsize)
        collect (list next-t x) into graph
        finally (return
                  (nconc graph (aif (ensure-last-step ode x (- next-t stepsize) t1)
                                 (list (list t1 it)))))))

(defun runge-kutta-values (ode x0 t0 t1 stepsize)
  "Calculate the values of the solution of the initial value problem
  at the points {t0, t0 + k * stepsize, t1}."
  (mapcar #'second (runge-kutta-graph ode x0 t0 t1 stepsize)))

(defun runge-kutta-values-standalone (ode x0 t0 t1 stepsize)
  "Calculate the values of the solution of the initial value problem
  at the points {t0, t0 + k * stepsize, t1}.
  This version does not calculate the whole graph first."
  (loop as next-t from t0 to t1 by stepsize
        for x = x0
        then (runge-kutta-single-step ode x next-t stepsize)
        collect x into vals
        finally (return (nconc vals (aif (ensure-last-step ode x (- next-t stepsize) t1)
                                      (list it))))))

(defun runge-kutta (ode x0 t0 t1 stepsize)
  "Calculate the value of the solution of the initial value problem
  at the point t1.
  This version does not calculate the whole graph first."
  (cadar (last (runge-kutta-graph ode x0 t0 t1 stepsize))))

(defun runge-kutta-standalone (ode x0 t0 t1 stepsize)
  "Calculate the value of the solution of the initial value problem
  at the point t1.
  This version does not calculate the whole graph first."
  (loop as next-t from t0 to t1 by stepsize
        for x = x0
        then (runge-kutta-single-step ode x next-t stepsize)
        finally (return (aif (ensure-last-step ode x (- next-t stepsize) t1)
                          it
                          x))))

(defun find-last-smaller (item alist)
  (loop for i = alist
        then (rest i)
        when (or (null (rest i))
                 (> (caar (rest i)) item))
        return i))

(defun cached-ode-solution (the-ode init-x0 init-t0 &optional (init-stepsize *stepsize*))
  "Create a closure caching already calculated values of the solution.
   Values between steps are interpolated on a straight line through adjecent known points.
   You can change the tableau between evolutions."
  (let ((stepsize (max *eps* init-stepsize))
        (known-graph (list (list init-t0 init-x0)))
        (furthest-step init-t0))
    (values
      (lambda (t1)
        (if (> t1 furthest-step)
          (let1 (graph (runge-kutta-graph the-ode (cadar (last known-graph)) furthest-step t1 stepsize))
            (appendf known-graph (rest graph))
            (setf furthest-step t1)
            (cadar (last graph)))
          (let1 (closest (find-last-smaller t1 known-graph))
            (destructuring-bind ((t0 x0) (t2 x2) &rest rst) closest
              (declare (ignore rst))
              (if (< t1 t0)
                (error "Can't go back in time."))
              (if (num-equal t0 t1 init-stepsize)
                x0
                (let1 (scale (/ (- t1 t0) (- t2 t0)))
                  (flet ((interpol (x0i x2i) (+ x0i (* scale x2i))))
                    (mapcar #'interpol x0 x2))))))))
      (lambda () known-graph))))
