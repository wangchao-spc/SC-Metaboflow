from SCMeTA.config import SYSTEM

RawFileReader = None
_raw_import_error = None

if SYSTEM == 'Darwin':
    try:
        from pyRawTools import MSLoader as RawFileReader
    except Exception as e:
        _raw_import_error = e
else:
    try:
        from .thermo_dotnet import MSLoader as RawFileReader
    except Exception as e:
        _raw_import_error = e

if RawFileReader is None:

    class _RawReaderUnavailable(object):
        """Placeholder so `import SCMeTA` never fails on machines where the
        Thermo .NET components could not be loaded (missing DLLs or no .NET
        Framework runtime). Raises a clear error only when .raw reading is
        actually attempted."""

        def __init__(self, *args, **kwargs):
            raise RuntimeError(
                "Thermo .RAW reading is unavailable on this computer: the "
                "ThermoFisher .NET components could not be loaded "
                f"({_raw_import_error}). This does not affect XLSX/CSV "
                "processing.")

    RawFileReader = _RawReaderUnavailable


__all__ = ["RawFileReader"]
