from typing import Sequence
import copy

import numpy as np


class Dot:
    def __init__(self,
                 x: float | int,
                 y: float | int,
                 color: str | None = None,
                 size: float | None = None,

                 ) -> None:
        self.x = x
        self.y = y
        self.color = color
        self.size = size


class Dots:
    """Stores dots with shared style properties (color, size)."""

    def __init__(self, color: str, size: float):
        self.x: np.ndarray = np.array([])
        self.y: np.ndarray = np.array([])
        self.color = color
        self.size = size

    def add_segment(self, x: float | int, y: float | int) -> None:
        """Add a line segment, inserting None to break between segments."""
        self.x = np.append(self.x, x)
        self.y = np.append(self.y, y)

    def join(self, other: "Dots") -> None:
        self.x = np.append(self.x, other.x)
        self.y = np.append(self.y, other.y)

    def matches(self, color: str, size: float | int) -> bool:
        """Return True if style matches."""
        return (
                self.color == color
                and self.size == size
        )


class Line:
    def __init__(self,
                 x: np.ndarray,
                 y: np.ndarray,
                 color: str | None = None,
                 width: float | None = None,
                 dash: str | None = None,
                 wave: bool = False
                 ) -> None:
        self.x = x
        self.y = y
        self.color = color
        self.width = width
        self.dash = dash
        self.wave = wave


class Lines:
    """Stores segments with shared style properties (color, width, dash, wave)."""

    def __init__(self, color: str, width: float, dash: str):
        self.x: np.ndarray = np.array([])
        self.y: np.ndarray = np.array([])
        self.color = color
        self.width = width
        self.dash = dash  # "solid", "dot", "dash", "longdash", "dashdot", or "longdashdot"

    def add_segment(self, x: np.ndarray, y: np.ndarray) -> None:
        """Add a line segment, inserting None to break between segments."""
        self.x = np.concatenate((self.x, x, [None]))
        self.y = np.concatenate((self.y, y, [None]))

    def join(self, other: "Lines") -> None:
        self.x = np.append(self.x, other.x)
        self.y = np.append(self.y, other.y)

    def matches(self, color: str, width: float, dash: str) -> bool:
        """Return True if style matches."""
        return (
                self.color == color
                and self.width == width
                and self.dash == dash
        )


class Fill:
    def __init__(self,
                 x: np.ndarray,
                 y: np.ndarray,
                 color: str
                 ) -> None:
        self.x = x
        self.y = y
        self.color = color


class Fills:
    """Stores filled polygon segments of the same color."""

    def __init__(self, color: str):
        self.x: np.ndarray = np.array([])
        self.y: np.ndarray = np.array([])
        self.color = color

    def add_segment(self, x: np.ndarray, y: np.ndarray) -> None:
        self.x = np.concatenate((self.x, x, [None]))
        self.y = np.concatenate((self.y, y, [None]))

    def join(self, other: "Fills") -> None:
        self.x = np.append(self.x, other.x)
        self.y = np.append(self.y, other.y)

    def matches(self, color: str) -> bool:
        return self.color == color


class DrawingContainer:
    """Unified container for lines, fills, and text, with Plotly export."""

    def __init__(self, layer_type: str | None = None):
        self.dots: list[Dots] = []
        self.lines: list[Lines] = []
        self.fills: list[Fills] = []

    def __str__(self):
        text = ""
        if self.dots:
            text += f"dots: {len(self.dots)} |"
        if self.lines:
            text += f"lines: {len(self.lines)} |"
        if self.fills:
            text += f"fills: {len(self.fills)} |"
        return text

    def is_empty(self) -> bool:
        return len(self.dots) == 0 and len(self.lines) == 0 and len(self.fills) == 0

    def add_dot(self, dot: Dot):
        for d in self.dots:
            if d.matches(dot.color, dot.size):
                d.add_segment(dot.x, dot.y)
                return

        new_dot = Dots(dot.color, dot.size)
        new_dot.add_segment(dot.x, dot.y)
        self.dots.append(new_dot)

    def add_line(self, line: Line):
        for l in self.lines:
            if l.matches(line.color, line.width, line.dash):
                l.add_segment(line.x, line.y)
                return

        new_line = Lines(line.color, line.width, line.dash)
        new_line.add_segment(line.x, line.y)
        self.lines.append(new_line)

    def add_fill(self, fill: Fill) -> None:
        for f in self.fills:
            if f.matches(fill.color):
                f.add_segment(fill.x, fill.y)
                return

        new_fill = Fills(fill.color)
        new_fill.add_segment(fill.x, fill.y)
        self.fills.append(new_fill)

    def add_objects(self, objs: Dot | Line | Fill | Sequence[Dot | Line | Fill ]):
        if not isinstance(objs, Sequence):
            objs = [objs]

        for obj in objs:
            if isinstance(obj, Line):
                self.add_line(obj)
            elif isinstance(obj, Fill):
                self.add_fill(obj)
            elif isinstance(obj, Dot):
                self.add_dot(obj)
            else:
                raise RuntimeError(f"Unknown type {type(obj)}")
