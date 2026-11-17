from PySide6.QtCore import QThread, Signal


class TCRAFTWorker(QThread):
    finished = Signal()
    error = Signal(str)
    log = Signal(str)

    def __init__(self, function, input_path, output_dir=None, **kwargs):
        super().__init__()

        self.function = function
        self.input_path = input_path
        self.output_dir = output_dir
        self.kwargs = kwargs

    def run(self):
        try:
            self.function(
                input_csv_path=self.input_path,
                output_dir=self.output_dir,
                **self.kwargs,
            )

        except Exception as exc:
            self.error.emit(str(exc))

        else:
            self.finished.emit()