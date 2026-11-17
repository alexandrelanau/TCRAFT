"""
Background execution of the TCRAFT engine for the GUI.

The engine reports its progress by printing to stdout. A GUI must never block
its event loop, so the run happens in a QThread and everything the engine prints
is forwarded to the window through Qt signals.
"""

import io
import re
import sys
import traceback

from PySide6.QtCore import QThread, Signal

#Intermediate tqdm progress bars would flood the log; only the completed one is
#worth keeping. Real progress is reported through the progress signal instead.
_TQDM_LINE = re.compile(r'^\s*(\d+)%\|')


class _SignalStream(io.TextIOBase):
    """
    A writable text stream that forwards whatever is written to it to a callable.
    Used to capture what the engine prints while it runs in the worker thread.
    """

    def __init__(self, emit):
        super().__init__()
        self._emit = emit
        self._buffer = ''

    def write(self, text):
        if not text:
            return 0

        self._buffer += text

        # Forward complete lines only, so the log stays readable. tqdm-style
        # carriage returns are treated as line breaks, so the FIRST of the two
        # delimiters ends the line; taking the last one would merge lines.
        while True:
            candidates = [i for i in (self._buffer.find('\n'), self._buffer.find('\r')) if i >= 0]
            if not candidates:
                break
            index = min(candidates)
            line = self._buffer[:index].rstrip('\r\n')
            self._buffer = self._buffer[index + 1:]

            self._emit_line(line)

        return len(text)

    def _emit_line(self, line):
        if not line.strip():
            return

        # Keep only the completed progress bar out of a tqdm series.
        match = _TQDM_LINE.match(line)
        if match and match.group(1) != '100':
            return

        self._emit(line)

    def flush(self):
        # tqdm calls flush() after every repaint of its progress bar, while the
        # line is still incomplete. Emitting here would leak every intermediate
        # bar into the log, so the pending text is kept until finish().
        pass

    def finish(self):
        """Emits whatever is left in the buffer at the end of the run."""
        if self._buffer:
            self._emit_line(self._buffer.rstrip('\r\n'))
        self._buffer = ''

    def isatty(self):
        return False


class TCRAFTWorker(QThread):
    """
    Runs one TCRAFT function in a background thread.

    Signals:
        log(str): one line printed by the engine.
        progress(int, int): TCRs processed, total TCRs.
        succeeded(str): the run finished; carries the output directory.
        failed(str, str): the run failed; carries the message and the traceback.
    """

    log = Signal(str)
    progress = Signal(int, int)
    succeeded = Signal(str)
    failed = Signal(str, str)

    def __init__(self, function, output_dir, supports_progress=False, **kwargs):
        super().__init__()

        self.function = function
        self.output_dir = output_dir
        self.supports_progress = supports_progress
        self.kwargs = kwargs

    def run(self):
        stream = _SignalStream(self.log.emit)

        # The engine writes with print(); capture it for the duration of the run.
        saved_stdout, saved_stderr = sys.stdout, sys.stderr
        sys.stdout = sys.stderr = stream

        try:
            kwargs = dict(self.kwargs)
            kwargs['output_dir'] = self.output_dir

            if self.supports_progress:
                kwargs['progress_callback'] = lambda done, total: self.progress.emit(done, total)

            self.function(**kwargs)

        except Exception as exc:
            stream.finish()
            self.failed.emit(str(exc) or exc.__class__.__name__, traceback.format_exc())

        else:
            stream.finish()
            self.succeeded.emit(self.output_dir)

        finally:
            sys.stdout, sys.stderr = saved_stdout, saved_stderr
