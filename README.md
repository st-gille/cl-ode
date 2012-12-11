Package CL-ODE
=====================

Solve systems of ordinary differential equations using Runge-Kutta methods.
Both explicit and implicit methods are possible, where implicit methdos are
solved with an implementation of Newtons method in C.

See `examples.lisp` for some examples.
Default Butcher-tableau is the explicit Euler method.
Wrap calls to runge-kutta\* in `with-tableau` to use a different method.
Choose a method via symbols from predefined tableaus
or define your own tableau with `make-butcher`.
Use `cached-ode-solution` for repeated evolutions.

This project uses [CFFI](https://github.com/cffi/cffi) to interface with a handcrafted C-library
and [Alexandria](http://common-lisp.net/project/alexandria/) for some general utilities.

Usage
-------------
Use `make CC_FLAGS=-DNDEBUG install` to install the shared library *libnewton*.
Compile a small test program calling the library with `make test`.
Load a sbcl repl with the package loaded with `make ode`.

Tested with `SBCL 1.0.55.0.debian` and `gcc version 4.6.3`.
