def db_setup(name=None):
    # Import the requested backend into a generic module object
    backend_name = 'backend_' + name
    backend_name = backend_name.lower()
    backend_name = 'nexthiv.db.%s' % backend_name
    # the last argument is specifies whether to use absolute or relative
    # imports. 0 means only perform absolute imports.
    backend_mod = __import__(backend_name, globals(), locals(), [backend_name], 0)
    return(backend_mod)
