import inspect
from types import SimpleNamespace


def workspace2struct_fun(exclude=None):
    frame = inspect.currentframe()
    if frame is None or frame.f_back is None:
        return SimpleNamespace()

    caller_locals = dict(frame.f_back.f_locals)
    if exclude is None:
        exclude = set()
    else:
        exclude = set(exclude)

    workspace = {
        name: value
        for name, value in caller_locals.items()
        if name not in exclude
    }
    return SimpleNamespace(**workspace)

