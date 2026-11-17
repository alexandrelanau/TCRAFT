"""
Entry point of the TCRAFT graphical interface.

Started by the TCRAFT-gui command, or directly by the bundled executable.
"""

import sys


def run_command_line():
    """
    The bundled executable is a GUI, but running it with arguments drives the
    engine from a terminal instead:

        TCRAFT generate tcrs.csv --output_dir results
        TCRAFT validate oligos.csv

    This also makes the built executable testable end to end.
    """
    command = sys.argv[1].lower()
    sys.argv = [f'TCRAFT-{command}'] + sys.argv[2:]

    if command == 'generate':
        from TCRAFT.Generate import run_CLI
    elif command == 'validate':
        from TCRAFT.Validate import run_CLI
    else:
        print(f'Unknown command "{command}". Expected "generate" or "validate".\n'
              f'Start TCRAFT without any argument to open the graphical interface.')
        sys.exit(2)

    run_CLI()


def main():
    # Arguments mean the user wants the command line, not the window.
    if len(sys.argv) > 1 and not sys.argv[1].startswith('-'):
        run_command_line()
        return

    from PySide6.QtWidgets import QApplication

    # Absolute import: PyInstaller runs this file as __main__, with no parent
    # package, so a relative import would fail in the bundled executable.
    from TCRAFT.gui.main_window import MainWindow

    app = QApplication(sys.argv)

    app.setApplicationName('TCRAFT')
    app.setApplicationDisplayName('TCRAFT')
    app.setOrganizationName('TCRAFT')

    window = MainWindow()
    window.show()

    sys.exit(app.exec())


if __name__ == '__main__':
    main()
