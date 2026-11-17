# -*- mode: python ; coding: utf-8 -*-
"""
PyInstaller build of the TCRAFT graphical interface.

Produces a single self-contained file (TCRAFT.exe on Windows, TCRAFT on Linux)
and, on macOS, a TCRAFT.app bundle. The reference data shipped in
TCRAFT/references is bundled with it, so the executable needs nothing else
installed on the machine.

Build with:   pyinstaller pyinstaller/TCRAFT.spec
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(SPECPATH).parent
PACKAGE_ROOT = PROJECT_ROOT / 'TCRAFT'
REFERENCES = PACKAGE_ROOT / 'references'

# The CSV/JSON/FASTA reference files the engine loads at import time. They are
# listed from the source tree rather than through collect_data_files(), which
# only finds them when TCRAFT is already installed in the build environment.
datas = [
    (str(path), 'TCRAFT/references')
    for path in sorted(REFERENCES.iterdir())
    if path.is_file() and path.suffix != '.pyc'
]

if not any(path.endswith('H_sapiens_codon_usage.json') for path, _ in datas):
    raise SystemExit(f'Reference data not found in {REFERENCES}; aborting the build.')

print(f'TCRAFT: bundling {len(datas)} reference files from {REFERENCES}')


a = Analysis(
    [str(PACKAGE_ROOT / 'gui' / 'main.py')],
    pathex=[str(PROJECT_ROOT)],
    binaries=[],
    datas=datas,
    # TCRAFT.references is only ever named as a string, in
    # importlib.resources.files('TCRAFT.references'), so the analysis cannot see
    # it. Without this the bundled app starts and then fails to find its data.
    hiddenimports=['TCRAFT.references'],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    # Nothing here is used by TCRAFT; leaving them out keeps the download smaller.
    excludes=[
        'tkinter',
        'matplotlib',
        'IPython',
        'jupyter',
        'notebook',
        'pytest',
        'PySide6.QtWebEngineCore',
        'PySide6.QtWebEngineWidgets',
        'PySide6.Qt3DCore',
        'PySide6.QtMultimedia',
        'PySide6.QtQuick',
        'PySide6.QtQml',
    ],
    noarchive=False,
)


pyz = PYZ(a.pure)


exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name='TCRAFT',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)


# On macOS the double-clickable artifact is an .app bundle, not a bare binary.
if sys.platform == 'darwin':
    app = BUNDLE(
        exe,
        name='TCRAFT.app',
        icon=None,
        bundle_identifier='edu.mit.birnbaumlab.tcraft',
        info_plist={
            'CFBundleName': 'TCRAFT',
            'CFBundleDisplayName': 'TCRAFT',
            'CFBundleShortVersionString': '1.2.0',
            'CFBundleVersion': '1.2.0',
            'NSHighResolutionCapable': True,
            # The app only ever reads and writes files the user picks itself.
            'LSApplicationCategoryType': 'public.app-category.education',
        },
    )
