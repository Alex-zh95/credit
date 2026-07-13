'''
Chart chrome and palette tokens: Kanagawa dark (rebelot/kanagawa.nvim).

The app commits to a single dark look. Sequential magnitude uses a one-hue
"wave" ramp built from Kanagawa's blue steps, anchored dark -> light so that
near-zero values recede into the dark chart surface.
'''

# Surfaces
PAGE = '#16161d'            # sumiInk0 — window plane
SURFACE = '#1f1f28'         # sumiInk1 — chart surface
SURFACE_RAISED = '#2a2a37'  # sumiInk2 — headers / controls
TOOLTIP_BG = '#223249'      # waveBlue1 — popup background

# Ink
INK = '#dcd7ba'             # fujiWhite — primary text
INK_SECONDARY = '#c8c093'   # oldWhite — secondary text (titles, status)
MUTED = '#727169'           # fujiGray — axis labels / ticks
GRID = '#363646'            # sumiInk3 — hairline gridlines
BASELINE = '#54546d'        # sumiInk4 — axis lines / borders

# Accents
ACCENT = '#7e9cd8'          # crystalBlue — hover marker
SELECTION = '#2d4f67'       # waveBlue2 — selection highlight
CRITICAL = '#ff5d62'        # peachRed — error status text (samuraiRed is
                            # sub-4.5:1 on sumiInk at status-bar text size)

# Sequential "wave" ramp, dark -> light (waveBlue1, waveBlue2, dragonBlue,
# springBlue, lightBlue). Monotone lightness; low values recede into SURFACE.
SEQ_WAVE = ['#223249', '#2d4f67', '#658594', '#7fb4ca', '#a3d4d5']


def sequential_cmap():
    '''Matplotlib colormap built from the Kanagawa wave ramp.'''
    from matplotlib.colors import LinearSegmentedColormap

    cmap = LinearSegmentedColormap.from_list('kanagawa_wave', SEQ_WAVE)
    cmap.set_bad(alpha=0.0)  # masked (no-data) cells render as gaps
    return cmap


def apply_qt_theme(app):
    '''
    Kanagawa the whole application: Fusion style + dark QPalette (this also
    makes matplotlib's toolbar invert its icons for a dark background) plus a
    small stylesheet for the pieces the palette does not reach.
    '''
    from PySide6.QtGui import QColor, QPalette

    app.setStyle('Fusion')
    p = QPalette()
    groups = (QPalette.ColorGroup.Active, QPalette.ColorGroup.Inactive)
    roles = {
        QPalette.ColorRole.Window: PAGE,
        QPalette.ColorRole.WindowText: INK,
        QPalette.ColorRole.Base: SURFACE,
        QPalette.ColorRole.AlternateBase: SURFACE_RAISED,
        QPalette.ColorRole.Text: INK,
        QPalette.ColorRole.Button: SURFACE_RAISED,
        QPalette.ColorRole.ButtonText: INK,
        QPalette.ColorRole.Highlight: SELECTION,
        QPalette.ColorRole.HighlightedText: INK,
        QPalette.ColorRole.ToolTipBase: TOOLTIP_BG,
        QPalette.ColorRole.ToolTipText: INK,
        QPalette.ColorRole.PlaceholderText: MUTED,
        QPalette.ColorRole.Light: BASELINE,
        QPalette.ColorRole.Midlight: GRID,
        QPalette.ColorRole.Mid: GRID,
        QPalette.ColorRole.Dark: PAGE,
        QPalette.ColorRole.Shadow: PAGE,
    }
    for group in groups:
        for role, colour in roles.items():
            p.setColor(group, role, QColor(colour))
    disabled = QPalette.ColorGroup.Disabled
    for role in (QPalette.ColorRole.WindowText, QPalette.ColorRole.Text,
                 QPalette.ColorRole.ButtonText):
        p.setColor(disabled, role, QColor(BASELINE))
    p.setColor(disabled, QPalette.ColorRole.Base, QColor(PAGE))
    app.setPalette(p)

    app.setStyleSheet(f'''
        QToolBar {{ background: {PAGE}; border: 0; spacing: 2px; }}
        QStatusBar {{ background: {PAGE}; color: {INK_SECONDARY}; }}
        QStatusBar::item {{ border: 0; }}
        QTableWidget {{
            background: {SURFACE};
            gridline-color: {GRID};
            border: 1px solid {GRID};
        }}
        QHeaderView::section {{
            background: {SURFACE_RAISED};
            color: {INK_SECONDARY};
            border: 0;
            padding: 3px 6px;
        }}
        QTableCornerButton::section {{ background: {SURFACE_RAISED}; border: 0; }}
    ''')
