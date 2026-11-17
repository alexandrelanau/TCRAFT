"""
TCRAFT main window.

Designed for a user who has never opened a terminal: pick a CSV, check what
TCRAFT understood from it, press one button, read the summary, open the folder.
"""

import os
from pathlib import Path

import pandas as pd

from PySide6.QtCore import Qt, QUrl
from PySide6.QtGui import QDesktopServices, QFont
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMainWindow,
    QMessageBox,
    QPlainTextEdit,
    QProgressBar,
    QPushButton,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from TCRAFT.Generate import (
    DEFAULT_SEED,
    REQUIRED_COLUMNS,
    generate_cdr3_oligos,
    summarize_input,
)
from TCRAFT.Validate import validate_cdr3_oligos

from .worker import TCRAFTWorker

AUTO_DETECT_LABEL = 'Auto-detect (recommended)'


class MainWindow(QMainWindow):

    def __init__(self):
        super().__init__()

        self.worker = None
        self.last_output_dir = None

        self.setWindowTitle('TCRAFT')
        self.resize(960, 760)

        self.build_ui()

    # ------------------------------------------------------------------ layout

    def build_ui(self):
        central = QWidget()
        layout = QVBoxLayout(central)

        title = QLabel('TCRAFT')
        title.setStyleSheet('font-size: 26px; font-weight: bold; padding: 6px 2px 0 2px;')

        subtitle = QLabel('TCR library oligo generation and assembly validation')
        subtitle.setStyleSheet('font-size: 13px; color: #666; padding: 0 2px 8px 2px;')

        layout.addWidget(title)
        layout.addWidget(subtitle)

        self.tabs = QTabWidget()
        self.tabs.addTab(self.create_generate_tab(), 'Generate oligos')
        self.tabs.addTab(self.create_validate_tab(), 'Validate assembly')
        layout.addWidget(self.tabs)

        self.progress = QProgressBar()
        self.progress.setTextVisible(True)
        self.progress.hide()
        layout.addWidget(self.progress)

        self.summary = QLabel()
        self.summary.setWordWrap(True)
        self.summary.setTextInteractionFlags(Qt.TextSelectableByMouse)
        self.summary.setStyleSheet(
            'background: #eef7ee; border: 1px solid #bcd9bc; border-radius: 4px; padding: 8px;'
        )
        self.summary.hide()
        layout.addWidget(self.summary)

        log_row = QHBoxLayout()
        log_label = QLabel('Log')
        log_label.setStyleSheet('font-weight: bold;')
        log_row.addWidget(log_label)
        log_row.addStretch()

        self.open_folder_button = QPushButton('Open results folder')
        self.open_folder_button.clicked.connect(self.open_results_folder)
        self.open_folder_button.setEnabled(False)
        log_row.addWidget(self.open_folder_button)

        layout.addLayout(log_row)

        self.log_output = QPlainTextEdit()
        self.log_output.setReadOnly(True)
        self.log_output.setMaximumBlockCount(5000)
        self.log_output.setPlaceholderText('TCRAFT output will appear here while it runs...')
        self.log_output.setFont(QFont('Monospace', 9))
        self.log_output.setMinimumHeight(180)
        layout.addWidget(self.log_output)

        self.setCentralWidget(central)

    def create_generate_tab(self):
        widget = QWidget()
        layout = QVBoxLayout(widget)

        input_group = QGroupBox('1. Input TCR list (CSV)')
        input_layout = QHBoxLayout(input_group)

        self.generate_input = QLineEdit()
        self.generate_input.setPlaceholderText('Select a CSV file containing your TCR sequences')
        self.generate_input.textChanged.connect(self.refresh_input_preview)

        input_button = QPushButton('Browse...')
        input_button.clicked.connect(lambda: self.select_file(self.generate_input))

        input_layout.addWidget(self.generate_input)
        input_layout.addWidget(input_button)
        layout.addWidget(input_group)

        self.preview = QLabel('Select an input CSV to see what TCRAFT reads from it.')
        self.preview.setWordWrap(True)
        self.preview.setTextInteractionFlags(Qt.TextSelectableByMouse)
        self.preview.setStyleSheet(
            'background: #f4f6f8; border: 1px solid #d7dde3; border-radius: 4px; padding: 8px; color: #333;'
        )
        layout.addWidget(self.preview)

        options_group = QGroupBox('2. Options')
        options_layout = QFormLayout(options_group)

        self.clonotype_combo = QComboBox()
        self.clonotype_combo.addItem(AUTO_DETECT_LABEL)
        self.clonotype_combo.setToolTip(
            'Column of your CSV holding the clonotype names kept next to each oligo.\n'
            'With Auto-detect, TCRAFT finds it on its own and generates TCR_0001... '
            'when your file has none.'
        )
        options_layout.addRow('Clonotype name column:', self.clonotype_combo)

        seed_row = QWidget()
        seed_layout = QHBoxLayout(seed_row)
        seed_layout.setContentsMargins(0, 0, 0, 0)

        self.seed_checkbox = QCheckBox('Reproducible run')
        self.seed_checkbox.setChecked(True)
        self.seed_checkbox.setToolTip(
            'With a fixed seed, the same input file always yields exactly the same oligos.'
        )
        self.seed_checkbox.toggled.connect(lambda on: self.seed_input.setEnabled(on))

        self.seed_input = QLineEdit(str(DEFAULT_SEED))
        self.seed_input.setMaximumWidth(110)
        self.seed_input.setToolTip('Change this value and run again if a few TCRs fail to generate.')

        seed_layout.addWidget(self.seed_checkbox)
        seed_layout.addWidget(QLabel('seed:'))
        seed_layout.addWidget(self.seed_input)
        seed_layout.addStretch()

        options_layout.addRow('Random seed:', seed_row)
        layout.addWidget(options_group)

        output_group = QGroupBox('3. Output directory')
        output_layout = QHBoxLayout(output_group)

        self.generate_output = QLineEdit()
        self.generate_output.setPlaceholderText('Filled in automatically once an input file is selected')

        output_button = QPushButton('Browse...')
        output_button.clicked.connect(lambda: self.select_directory(self.generate_output))

        output_layout.addWidget(self.generate_output)
        output_layout.addWidget(output_button)
        layout.addWidget(output_group)

        self.generate_button = QPushButton('Generate CDR3 oligos')
        self.generate_button.setMinimumHeight(46)
        self.generate_button.setStyleSheet('font-size: 15px; font-weight: bold;')
        self.generate_button.clicked.connect(self.run_generate)
        layout.addWidget(self.generate_button)

        layout.addStretch()
        return widget

    def create_validate_tab(self):
        widget = QWidget()
        layout = QVBoxLayout(widget)

        explanation = QLabel(
            'Checks that oligos produced by TCRAFT are error-free and assemble, by Golden Gate, '
            'into the intended TCR sequences. Give it one of the oligo CSV files from a Generate run.'
        )
        explanation.setWordWrap(True)
        explanation.setStyleSheet('color: #555; padding-bottom: 4px;')
        layout.addWidget(explanation)

        input_group = QGroupBox('1. Oligo file produced by TCRAFT (CSV)')
        input_layout = QHBoxLayout(input_group)

        self.validate_input = QLineEdit()
        self.validate_input.setPlaceholderText('Select an oligo CSV generated by TCRAFT')
        self.validate_input.textChanged.connect(self.refresh_validate_output)

        input_button = QPushButton('Browse...')
        input_button.clicked.connect(lambda: self.select_file(self.validate_input))

        input_layout.addWidget(self.validate_input)
        input_layout.addWidget(input_button)
        layout.addWidget(input_group)

        output_group = QGroupBox('2. Output directory')
        output_layout = QHBoxLayout(output_group)

        self.validate_output = QLineEdit()
        self.validate_output.setPlaceholderText('Filled in automatically once an input file is selected')

        output_button = QPushButton('Browse...')
        output_button.clicked.connect(lambda: self.select_directory(self.validate_output))

        output_layout.addWidget(self.validate_output)
        output_layout.addWidget(output_button)
        layout.addWidget(output_group)

        self.validate_button = QPushButton('Validate assembly')
        self.validate_button.setMinimumHeight(46)
        self.validate_button.setStyleSheet('font-size: 15px; font-weight: bold;')
        self.validate_button.clicked.connect(self.run_validate)
        layout.addWidget(self.validate_button)

        layout.addStretch()
        return widget

    # ------------------------------------------------------------------ inputs

    def select_file(self, line_edit):
        filename, _ = QFileDialog.getOpenFileName(
            self, 'Select CSV file', line_edit.text() or '', 'CSV files (*.csv);;All files (*)'
        )
        if filename:
            line_edit.setText(filename)

    def select_directory(self, line_edit):
        directory = QFileDialog.getExistingDirectory(
            self, 'Select output directory', line_edit.text() or ''
        )
        if directory:
            line_edit.setText(directory)

    def default_output_dir(self, input_path, prefix):
        """
        Puts results next to the input file, so a user who never chooses a
        directory still knows where their files are.
        """
        date_str = pd.Timestamp.now().strftime('%Y%m%d')
        return str(Path(input_path).resolve().parent / f'TCRAFT-{prefix}_{date_str}')

    def refresh_validate_output(self):
        path = self.validate_input.text().strip()
        if path and Path(path).is_file():
            self.validate_output.setText(self.default_output_dir(path, 'validate'))

    def refresh_input_preview(self):
        """
        Reads the selected CSV and reports what TCRAFT understood from it, before
        anything is generated.
        """
        path = self.generate_input.text().strip()

        self.clonotype_combo.blockSignals(True)
        self.clonotype_combo.clear()
        self.clonotype_combo.addItem(AUTO_DETECT_LABEL)
        self.clonotype_combo.blockSignals(False)

        if not path:
            self.preview.setText('Select an input CSV to see what TCRAFT reads from it.')
            return

        if not Path(path).is_file():
            self.preview.setText('<b>This file does not exist.</b>')
            return

        self.generate_output.setText(self.default_output_dir(path, 'generate'))

        info = summarize_input(path)

        if info['error']:
            self.preview.setText(f"<b>{info['error']}</b>")
            return

        for column in info['columns']:
            self.clonotype_combo.addItem(column)

        if info['clonotype_column'] is not None:
            index = self.clonotype_combo.findText(info['clonotype_column'])
            if index > 0:
                self.clonotype_combo.setItemText(
                    0, f"{AUTO_DETECT_LABEL} - \"{info['clonotype_column']}\""
                )

        lines = [f"<b>{info['n_rows']}</b> rows read, <b>{len(info['columns'])}</b> columns."]

        if info['missing_columns']:
            lines.append(
                '<b style="color:#a33;">Missing required column(s): '
                + ', '.join(info['missing_columns'])
                + '</b><br>Required: '
                + ', '.join(REQUIRED_COLUMNS)
            )
        else:
            lines.append('All six required TCR columns are present.')

            if info['clonotype_column'] is not None:
                lines.append(
                    f"Clonotype names will be taken from column <b>\"{info['clonotype_column']}\"</b>."
                )
            else:
                lines.append(
                    'No clonotype name column found: names <b>TCR_0001...</b> will be generated, '
                    'so every oligo still carries a name.'
                )

            if info['n_duplicates']:
                lines.append(
                    f"<b>{info['n_duplicates']}</b> duplicate TCR(s) found and will be removed "
                    f"(listed in All_TCR_removed.csv); <b>{info['n_unique']}</b> oligos will be generated."
                )
            else:
                lines.append('No duplicate TCR found.')

        self.preview.setText('<br>'.join(lines))

    # --------------------------------------------------------------------- run

    def read_seed(self):
        if not self.seed_checkbox.isChecked():
            return None

        text = self.seed_input.text().strip()

        if not text:
            return DEFAULT_SEED

        try:
            return int(text)
        except ValueError:
            raise ValueError(f'The random seed must be a whole number, not "{text}".')

    def run_generate(self):
        input_path = self.generate_input.text().strip()

        if not input_path:
            QMessageBox.warning(self, 'Missing input', 'Please select a CSV input file.')
            return

        if not Path(input_path).is_file():
            QMessageBox.warning(self, 'Invalid input', 'The selected input file does not exist.')
            return

        info = summarize_input(input_path)

        if info['error']:
            QMessageBox.critical(self, 'Invalid input', info['error'])
            return

        if info['missing_columns']:
            QMessageBox.critical(
                self,
                'Missing columns',
                'This CSV is missing the required column(s):\n\n'
                + ', '.join(info['missing_columns'])
                + '\n\nRequired columns are:\n'
                + ', '.join(REQUIRED_COLUMNS),
            )
            return

        try:
            seed = self.read_seed()
        except ValueError as exc:
            QMessageBox.warning(self, 'Invalid seed', str(exc))
            return

        clonotype_column = None
        if self.clonotype_combo.currentIndex() > 0:
            clonotype_column = self.clonotype_combo.currentText()

        output_dir = self.generate_output.text().strip() or self.default_output_dir(
            input_path, 'generate'
        )

        self.start_worker(
            generate_cdr3_oligos,
            input_path,
            output_dir,
            supports_progress=True,
            total=info['n_unique'] or info['n_rows'],
            input_csv_path=input_path,
            clonotype_id_column=clonotype_column,
            seed=seed,
        )

    def run_validate(self):
        input_path = self.validate_input.text().strip()

        if not input_path:
            QMessageBox.warning(self, 'Missing input', 'Please select a CSV input file.')
            return

        if not Path(input_path).is_file():
            QMessageBox.warning(self, 'Invalid input', 'The selected input file does not exist.')
            return

        output_dir = self.validate_output.text().strip() or self.default_output_dir(
            input_path, 'validate'
        )

        self.start_worker(
            validate_cdr3_oligos,
            input_path,
            output_dir,
            supports_progress=False,
            total=0,
            input_csv_path=input_path,
        )

    def start_worker(self, function, input_path, output_dir, supports_progress, total, **kwargs):
        self.log_output.clear()
        self.summary.hide()
        self.open_folder_button.setEnabled(False)
        self.last_output_dir = output_dir

        self.log_output.appendPlainText(f'Input:  {input_path}')
        self.log_output.appendPlainText(f'Output: {output_dir}')
        self.log_output.appendPlainText('')

        if total:
            self.progress.setRange(0, total)
            self.progress.setValue(0)
            self.progress.setFormat('%v / %m TCRs')
        else:
            self.progress.setRange(0, 0)
            self.progress.setFormat('Working...')

        self.set_running(True)

        self.worker = TCRAFTWorker(
            function, output_dir, supports_progress=supports_progress, **kwargs
        )

        self.worker.log.connect(self.append_log)
        self.worker.progress.connect(self.on_progress)
        self.worker.succeeded.connect(self.on_worker_succeeded)
        self.worker.failed.connect(self.on_worker_failed)
        self.worker.start()

    def set_running(self, running):
        self.generate_button.setEnabled(not running)
        self.validate_button.setEnabled(not running)
        self.tabs.setEnabled(not running)

        self.generate_button.setText('Generating, please wait...' if running else 'Generate CDR3 oligos')

        if running:
            self.progress.show()
        else:
            self.progress.hide()

    def append_log(self, line):
        self.log_output.appendPlainText(line)

    def on_progress(self, done, total):
        if self.progress.maximum() != total:
            self.progress.setRange(0, total)
        self.progress.setValue(done)

    # ----------------------------------------------------------------- results

    def read_metadata(self, output_dir):
        """
        Turns Oligo_Pool_Metadata.csv into the lines shown in the summary box.
        """
        path = os.path.join(output_dir, 'Oligo_Pool_Metadata.csv')

        if not os.path.exists(path):
            return []

        try:
            table = pd.read_csv(path)
        except Exception:
            return []

        lines = []
        for _, row in table.iterrows():
            name, value = str(row.iloc[0]), str(row.iloc[1])
            if value in ('nan', ''):
                continue
            if name.startswith('Avg Oligo Size'):
                try:
                    value = f'{float(value):.0f} bp'
                except ValueError:
                    pass
            lines.append(f'{name}: {value}')

        return lines

    def on_worker_succeeded(self, output_dir):
        self.set_running(False)
        self.open_folder_button.setEnabled(True)

        self.log_output.appendPlainText('')
        self.log_output.appendPlainText('TCRAFT completed successfully.')

        lines = ['<b>Done.</b> Results saved in:', f'<code>{output_dir}</code>']

        metadata = self.read_metadata(output_dir)
        if metadata:
            lines.append('')
            lines.extend(metadata)

        report = os.path.join(output_dir, 'Padding_report.csv')
        if os.path.exists(report):
            try:
                if not pd.read_csv(report).empty:
                    lines.append('')
                    lines.append('Some oligos needed special padding: see Padding_report.csv.')
            except Exception:
                pass

        self.summary.setText('<br>'.join(lines))
        self.summary.setStyleSheet(
            'background: #eef7ee; border: 1px solid #bcd9bc; border-radius: 4px; padding: 8px;'
        )
        self.summary.show()

    def on_worker_failed(self, message, trace):
        self.set_running(False)

        self.log_output.appendPlainText('')
        self.log_output.appendPlainText(f'ERROR: {message}')
        self.log_output.appendPlainText('')
        self.log_output.appendPlainText(trace)

        if self.last_output_dir and os.path.isdir(self.last_output_dir):
            self.open_folder_button.setEnabled(True)

        self.summary.setText(f'<b>The run failed.</b><br>{message}<br><br>See the log below for details.')
        self.summary.setStyleSheet(
            'background: #fdeeee; border: 1px solid #e0b4b4; border-radius: 4px; padding: 8px;'
        )
        self.summary.show()

        QMessageBox.critical(self, 'TCRAFT error', message)

    def open_results_folder(self):
        if self.last_output_dir and os.path.isdir(self.last_output_dir):
            QDesktopServices.openUrl(QUrl.fromLocalFile(self.last_output_dir))
        else:
            QMessageBox.warning(self, 'TCRAFT', 'The results folder does not exist (yet).')
