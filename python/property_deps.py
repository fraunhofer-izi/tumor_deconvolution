def del_if_there(instance, varnames):
    """Deletes all variables in varnames if they exist in instance."""
    for var in varnames:
        try:
            delattr(instance, var)
        except AttributeError:
            pass

def property_deps(*args, verbose=False, setter_hook=None):
    """

    A function decorator to generate set, get and delete methods for a
    property on which others depend.

    Dependent variables are deleted when the property changes or is deleted.
    The decorated function will be used as getter that is executed once at
    request if there is no stored value and cached until it is deleted or
    set manually.
 
    Name and docstring of the decoarted function are preserved in the property.

    Parameters
    ----------
    args : string
        Each unnamed argument must be a variable names of a depending variable.
    verbose : bool (optional)
        Print message when the getter function is called.
    setter_hook : callalble (optional)
        A function that is called before the new value is set.
        It recives the instance and the new value as argument:
        `setter_hook(self, value)`

    Returns
    -------
    dep_property : class
        A class with predifined __get__, __set__ and __delete__ methods.

    Usage
    -----
    class some_class(some_super_class):

        def _my_hook(self, value):
            print(type(self)) # do somethin with self

        @property_deps('some', 'dependent', 'variables', setter_hook=_my_hook)
        def my_property(self):
            ...
            return value
    """

    assert isinstance(verbose, bool), '`verbose` musst be boolean.'
    if setter_hook is None:
        def setter_hook(*args):
            pass
    assert callable(setter_hook), '`setter_hook` musst be callable.'

    class baseClass:
        def __init__(self, func):
            self.getter = func
            self.fname = '_cashed_dep_prop_' + func.__name__
        def __get__(self, instance, owner):
            if instance is None:
                return self
            if self.getter is None:
                raise AttributeError("unreadable attribute")
            if not hasattr(instance, self.fname):
                if verbose:
                    print('Getting {} ...'.format(self.getter.__name__))
                setattr(instance, self.fname, self.getter(instance))
            return getattr(instance, self.fname)
        def __set__(self, instance, value):
            setter_hook(instance, value)
            del_if_there(instance, args)
            setattr(instance, self.fname, value)
        def __delete__(self, instance):
            del_if_there(instance, (self.fname,) + args)

    def dep_property(func):
        clsdict = {'__doc__': func.__doc__}
        theClass = type(func.__name__, (baseClass,), clsdict)
        return theClass(func)

    return dep_property
