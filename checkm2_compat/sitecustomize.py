"""sitecustomize injected for CheckM2 1.1.0 + Python 3.12 compatibility.

CheckM2 launches `mp.Process(target=self.__set_up_prodigal_thread, ...)`
and similar. The bound method `self.__set_up_prodigal_thread` resolves to
`self._Predictor__set_up_prodigal_thread` after Python name mangling, so
that name is what gets pickled onto the process's run queue.

On Python 3.12 the multiprocessing pickling protocol changed: when the
worker reconstructs the bound method it asks for the *unmangled* name
`__set_up_prodigal_thread`, which doesn't exist on the class. Result:

    AttributeError: 'Predictor' object has no attribute
    '__set_up_prodigal_thread'. Did you mean:
    '_Predictor__set_up_prodigal_thread'?

This sitecustomize registers an audit hook so that, the moment
`checkm2.predictQuality` finishes loading, we add unmangled aliases for
the four mangled Predictor methods that CheckM2 uses as Process targets.
It also forces `multiprocessing.Pool` to use the `fork` start method, so
when CheckM2 falls back to a Pool (e.g. for the DIAMOND step) the
private-method pickling is bypassed entirely.
"""
import multiprocessing as _mp
import sys


def _install_fork_pool_patch() -> None:
    if getattr(_mp.Pool, "_db_builder_cli_fork_patch", False):
        return
    _orig_pool = _mp.Pool

    def _wrapper(*args, **kwargs):
        if "context" not in kwargs:
            kwargs["context"] = _mp.get_context("fork")
        return _orig_pool(*args, **kwargs)

    _wrapper._db_builder_cli_fork_patch = True  # type: ignore[attr-defined]
    _mp.Pool = _wrapper


def _patch_predictor_class() -> None:
    module = sys.modules.get("checkm2.predictQuality")
    if module is None:
        return
    P = getattr(module, "Predictor", None)
    if P is None:
        return
    for m in (
        "__set_up_prodigal_thread",
        "__set_up_metadata_thread",
        "__reportProgress",
        "__report_progress_metadata",
    ):
        mangled = "_Predictor" + m
        if hasattr(P, mangled) and not hasattr(P, m):
            setattr(P, m, getattr(P, mangled))


def _audit(event: str, args: tuple) -> None:
    if event == "import" and args and args[0] == "checkm2.predictQuality":
        # Defer one event loop turn so the module finishes loading.
        try:
            _patch_predictor_class()
        except Exception:
            pass


_install_fork_pool_patch()
sys.addaudithook(_audit)
