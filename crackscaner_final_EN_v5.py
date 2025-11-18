#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
svg_viewer_boundary_scale_axes_choice_axisflip_scanline_multi_angle.py

拡張版:
- これまでの機能（水平単一走査線・水平複数走査線・スケール・ミニ座標軸・Convergence plot別ウインドウ・Export CSV）はそのまま。
- 追加機能:
    * 0〜180°未満を一定Angle step（15° or 10°）で回転した走査線解析
    * 各角度ごとに「複数走査線」で指定した本数の走査線を等間隔で発生
    * 各角度について線分図＋交点を描画（メイン画面のタブで角度切り替え）
    * Convergence plotは別ウインドウに表示し、メイン画面でタブを切り替えると
      その角度の収束曲線にUpdate
    * 全角度・全走査線の結果を 1 つの CSV に出力（角度列を追加）
    * Rose diagram（0〜360°表示／シンメトリー考慮）を別ウインドウに表示
      - 横軸の「総走査線長さ L_total」を指定して、その時点の平均密度を各角度で評価
      @2025 Takato Takemura, Geomechanics Lab., Nihon Univ.
"""

import sys
import os
import re
import math
import csv
import tkinter as tk
from tkinter import filedialog, messagebox
import tkinter.ttk as ttk
from xml.etree import ElementTree as ET

import numpy as np
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

# --- 日本語が文字化けしないようにフォントを指定 ---
matplotlib.rcParams["font.family"] = [
    "Hiragino Sans",      # Macの標準日本語フォント
    "Yu Gothic",          # Win, Mac両方にあることが多い
    "Meiryo",             # Windows
    "MS Gothic"           # Windows古い環境用
]
matplotlib.rcParams["axes.unicode_minus"] = False  # マイナス記号を正常表示


# ---------------------- Parse helpers ----------------------
def _parse_points_attr(points_str: str):
    s = re.sub(r"[,\s]+", " ", points_str.strip())
    vals = [float(v) for v in s.split() if v.strip()]
    if len(vals) % 2 != 0:
        vals = vals[:-1]
    return list(zip(vals[0::2], vals[1::2]))


def _cubic_bezier(p0, p1, p2, p3, t):
    x = ((1 - t) ** 3) * p0[0] + 3 * ((1 - t) ** 2) * t * p1[0] + 3 * (1 - t) * (t ** 2) * p2[0] + (t ** 3) * p3[0]
    y = ((1 - t) ** 3) * p0[1] + 3 * ((1 - t) ** 2) * t * p1[1] + 3 * (1 - t) * (t ** 2) * p2[1] + (t ** 3) * p3[1]
    return (x, y)


def _quadratic_bezier(p0, p1, p2, t):
    x = (1 - t) ** 2 * p0[0] + 2 * (1 - t) * t * p1[0] + t ** 2 * p2[0]
    y = (1 - t) ** 2 * p0[1] + 2 * (1 - t) * t * p1[1] + t ** 2 * p2[1]
    return (x, y)


def _tokenize_path_d(d: str):
    # ★ S / s もコマンドとして扱うように修正
    return re.findall(r"[MmLlHhVvCcQqSsZz]|[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", d)


def _path_to_polyline(d: str, steps_curve: int = 32):
    tokens = _tokenize_path_d(d)
    i = 0
    x = y = 0.0
    start_x = start_y = 0.0
    last_cmd = None
    # ★ 直前の cubic Bézier の第2制御点（C / c / S / s 用）
    last_cx = last_cy = None

    polylines = []
    current_poly = []

    def flush_poly():
        nonlocal current_poly
        if len(current_poly) >= 2:
            polylines.append(current_poly)
        current_poly = []

    def add_point(px, py):
        nonlocal current_poly
        if (not current_poly) or (current_poly[-1][0] != px or current_poly[-1][1] != py):
            current_poly.append((px, py))

    while i < len(tokens):
        tok = tokens[i]
        i += 1
        if re.match(r"[MmLlHhVvCcQqSsZz]", tok):
            cmd = tok
        else:
            if last_cmd is None:
                break
            cmd = last_cmd
            i -= 1

        # --- M/m ---
        if cmd in "Mm":
            rel = (cmd == "m")
            if i + 1 >= len(tokens):
                break
            x0 = float(tokens[i])
            y0 = float(tokens[i + 1])
            i += 2
            if rel:
                x += x0
                y += y0
            else:
                x, y = x0, y0
            flush_poly()
            add_point(x, y)
            start_x, start_y = x, y
            last_cx = last_cy = None
            while i + 1 < len(tokens) and not re.match(r"[MmLlHhVvCcQqSsZz]", tokens[i]):
                x1 = float(tokens[i])
                y1 = float(tokens[i + 1])
                i += 2
                if rel:
                    x += x1
                    y += y1
                else:
                    x, y = x1, y1
                add_point(x, y)

        # --- L/l ---
        elif cmd in "Ll":
            rel = (cmd == "l")
            while i + 1 < len(tokens) and not re.match(r"[MmLlHhVvCcQqSsZz]", tokens[i]):
                dx = float(tokens[i])
                dy = float(tokens[i + 1])
                i += 2
                if rel:
                    x += dx
                    y += dy
                else:
                    x, y = dx, dy
                add_point(x, y)
            last_cx = last_cy = None

        # --- H/h ---
        elif cmd in "Hh":
            rel = (cmd == "h")
            while i < len(tokens) and not re.match(r"[MmLlHhVvCcQqSsZz]", tokens[i]):
                vx = float(tokens[i])
                i += 1
                x = x + vx if rel else vx
                add_point(x, y)
            last_cx = last_cy = None

        # --- V/v ---
        elif cmd in "Vv":
            rel = (cmd == "v")
            while i < len(tokens) and not re.match(r"[MmLlHhVvCcQqSsZz]", tokens[i]):
                vy = float(tokens[i])
                i += 1
                y = y + vy if rel else vy
                add_point(x, y)
            last_cx = last_cy = None

        # --- C/c (cubic Bézier) ---
        elif cmd in "Cc":
            rel = (cmd == "c")
            while i + 5 < len(tokens) and not re.match(r"[MmLlHhVvCcQqSsZz]", tokens[i]):
                x1 = float(tokens[i])
                y1 = float(tokens[i + 1])
                x2 = float(tokens[i + 2])
                y2 = float(tokens[i + 3])
                x3 = float(tokens[i + 4])
                y3 = float(tokens[i + 5])
                i += 6

                p0 = (x, y)
                if rel:
                    p1 = (x + x1, y + y1)
                    p2 = (x + x2, y + y2)
                    p3 = (x + x3, y + y3)
                    x, y = p3
                else:
                    p1 = (x1, y1)
                    p2 = (x2, y2)
                    p3 = (x3, y3)
                    x, y = p3

                # 直前の control point(第2制御点) を保存 → S/s で使用
                last_cx, last_cy = p2

                for k in range(1, steps_curve + 1):
                    t = k / steps_curve
                    add_point(*_cubic_bezier(p0, p1, p2, p3, t))

        # --- S/s (smooth cubic Bézier) ★追加 ---
        elif cmd in "Ss":
            rel = (cmd == "s")
            while i + 3 < len(tokens) and not re.match(r"[MmLlHhVvCcQqSsZz]", tokens[i]):
                x2 = float(tokens[i])
                y2 = float(tokens[i + 1])
                x3 = float(tokens[i + 2])
                y3 = float(tokens[i + 3])
                i += 4

                p0 = (x, y)

                # 直前が C/c or S/s のときだけ制御点を反転
                if last_cmd in "CcSs" and last_cx is not None:
                    p1x = 2 * x - last_cx
                    p1y = 2 * y - last_cy
                else:
                    p1x, p1y = x, y

                if rel:
                    p1 = (p1x, p1y)
                    p2 = (x + x2, y + y2)
                    p3 = (x + x3, y + y3)
                    x, y = p3
                else:
                    p1 = (p1x, p1y)
                    p2 = (x2, y2)
                    p3 = (x3, y3)
                    x, y = p3

                last_cx, last_cy = p2

                for k in range(1, steps_curve + 1):
                    t = k / steps_curve
                    add_point(*_cubic_bezier(p0, p1, p2, p3, t))

        # --- Q/q (quadratic Bézier) ---
        elif cmd in "Qq":
            rel = (cmd == "q")
            while i + 3 < len(tokens) and not re.match(r"[MmLlHhVvCcQqSsZz]", tokens[i]):
                x1 = float(tokens[i])
                y1 = float(tokens[i + 1])
                x2 = float(tokens[i + 2])
                y2 = float(tokens[i + 3])
                i += 4

                p0 = (x, y)
                if rel:
                    p1 = (x + x1, y + y1)
                    p2 = (x + x2, y + y2)
                    x, y = p2
                else:
                    p1 = (x1, y1)
                    p2 = (x2, y2)
                    x, y = p2

                last_cx = last_cy = None

                for k in range(1, steps_curve + 1):
                    t = k / steps_curve
                    add_point(*_quadratic_bezier(p0, p1, p2, t))

        # --- Z/z ---
        elif cmd in "Zz":
            add_point(start_x, start_y)
            if current_poly:
                polylines.append(current_poly)
            current_poly = []
            x, y = start_x, start_y
            last_cx = last_cy = None

        last_cmd = cmd

    if current_poly:
        polylines.append(current_poly)
    return polylines, ('Z' in ''.join(tokens) or 'z' in ''.join(tokens))


# ---------------------- Extraction ----------------------
def extract_all(svg_text: str):
    svg_text_clean = re.sub(r'xmlns="[^"]+"', '', svg_text, count=1)
    root = ET.fromstring(svg_text_clean)

    line_segments = []
    boundary_segments = []
    bxs, bys = [], []

    for elem in root.iter():
        tag = elem.tag.split('}')[-1]

        if tag == "line":
            try:
                x1 = float(elem.attrib.get("x1", "nan"))
                y1 = float(elem.attrib.get("y1", "nan"))
                x2 = float(elem.attrib.get("x2", "nan"))
                y2 = float(elem.attrib.get("y2", "nan"))
                if all(math.isfinite(v) for v in (x1, y1, x2, y2)):
                    line_segments.append((x1, y1, x2, y2))
            except Exception:
                pass

        elif tag == "polyline":
            pts_str = elem.attrib.get("points", "").strip()
            if pts_str:
                try:
                    pts = _parse_points_attr(pts_str)
                    for (x1, y1), (x2, y2) in zip(pts, pts[1:]):
                        if all(math.isfinite(v) for v in (x1, y1, x2, y2)):
                            line_segments.append((x1, y1, x2, y2))
                except Exception:
                    pass

        elif tag == "polygon":
            pts_str = elem.attrib.get("points", "").strip()
            if pts_str:
                try:
                    pts = _parse_points_attr(pts_str)
                    if pts and pts[0] != pts[-1]:
                        pts.append(pts[0])
                    for (x1, y1), (x2, y2) in zip(pts, pts[1:]):
                        if all(math.isfinite(v) for v in (x1, y1, x2, y2)):
                            boundary_segments.append((x1, y1, x2, y2))
                            bxs.extend([x1, x2])
                            bys.extend([y1, y2])
                except Exception:
                    pass

        elif tag == "rect":
            try:
                x = float(elem.attrib.get("x", "0"))
                y = float(elem.attrib.get("y", "0"))
                w = float(elem.attrib.get("width", "0"))
                h = float(elem.attrib.get("height", "0"))
                p = [(x, y), (x + w, y), (x + w, y + h), (x, y + h), (x, y)]
                for (x1, y1), (x2, y2) in zip(p, p[1:]):
                    if all(math.isfinite(v) for v in (x1, y1, x2, y2)):
                        boundary_segments.append((x1, y1, x2, y2))
                        bxs.extend([x1, x2])
                        bys.extend([y1, y2])
            except Exception:
                pass

        elif tag == "path":
            d = elem.attrib.get("d", "").strip()
            if d:
                try:
                    polylines, closed = _path_to_polyline(d, steps_curve=32)
                    if closed:
                        for poly in polylines:
                            if len(poly) >= 2 and poly[0] != poly[-1]:
                                poly = poly + [poly[0]]
                            for (x1, y1), (x2, y2) in zip(poly, poly[1:]):
                                if all(math.isfinite(v) for v in (x1, y1, x2, y2)):
                                    boundary_segments.append((x1, y1, x2, y2))
                                    bxs.extend([x1, x2])
                                    bys.extend([y1, y2])
                    else:
                        for poly in polylines:
                            for (x1, y1), (x2, y2) in zip(poly, poly[1:]):
                                if all(math.isfinite(v) for v in (x1, y1, x2, y2)):
                                    line_segments.append((x1, y1, x2, y2))
                except Exception:
                    # パス解析で例外が出た場合は、そのパスは無視
                    pass

    boundary_bbox = None    # (xmin, xmax, ymin, ymax)
    if bxs and bys:
        boundary_bbox = (min(bxs), max(bxs), min(bys), max(bys))
    return line_segments, boundary_segments, boundary_bbox



def compute_crack_lengths(svg_text: str):
    """
    SVG 内の <line> および <polyline> 要素を「クラック」とみなし，
    1 要素 = 1 本として端点間距離を長さとして返すヘルパー関数。
    - <line>: (x1, y1)–(x2, y2) の距離
    - <polyline>: 最初と最後の点の距離
    他の要素（polygon, rect, path など）はクラックには含めない。
    """
    svg_text_clean = re.sub(r'xmlns="[^"]+"', '', svg_text, count=1)
    root = ET.fromstring(svg_text_clean)

    crack_lengths = []

    for elem in root.iter():
        tag = elem.tag.split('}')[-1]

        if tag == "line":
            try:
                x1 = float(elem.attrib.get("x1", "nan"))
                y1 = float(elem.attrib.get("y1", "nan"))
                x2 = float(elem.attrib.get("x2", "nan"))
                y2 = float(elem.attrib.get("y2", "nan"))
                if all(math.isfinite(v) for v in (x1, y1, x2, y2)):
                    L = math.hypot(x2 - x1, y2 - y1)
                    crack_lengths.append(L)
            except Exception:
                pass

        elif tag == "polyline":
            pts_str = elem.attrib.get("points", "").strip()
            if pts_str:
                try:
                    pts = _parse_points_attr(pts_str)
                    if len(pts) >= 2:
                        (x1, y1) = pts[0]
                        (x2, y2) = pts[-1]
                        if all(math.isfinite(v) for v in (x1, y1, x2, y2)):
                            L = math.hypot(x2 - x1, y2 - y1)
                            crack_lengths.append(L)
                except Exception:
                    pass

    return crack_lengths


# ---------------------- Geometry helpers ----------------------
def _intersect_segment_with_horizontal(x1, y1, x2, y2, y0, eps=1e-9):
    """
    水平線 y=y0 用（従来の関数）
    """
    if abs(y2 - y1) < eps:
        return None
    t = (y0 - y1) / (y2 - y1)
    if t < -eps or t > 1 + eps:
        return None
    t = max(0.0, min(1.0, t))
    x = x1 + t * (x2 - x1)
    return x


def analyze_scanline_general(line_segments, boundary_segments, d, n, c, eps=1e-9):
    """
    一般角度の走査線解析
    line_segments, boundary_segments: (x1,y1,x2,y2) のリスト
    d: 走査線方向の単位ベクトル (dx,dy)
    n: その法線 (-sinθ, cosθ) など
    c: n・p = c で定義される直線
    戻り値: dict または None
    dict:
        {
          'pos': c,
          'length': L,
          'hits': N,
          'density': N/L,
          'boundary_points': [ (p_min), (p_max) ],
          'hit_points': [ (xi, yi), ... ]
        }
    """
    # 境界との交点
    b_points = []
    for (x1, y1, x2, y2) in boundary_segments:
        p1 = (x1, y1)
        p2 = (x2, y2)
        vx = p2[0] - p1[0]
        vy = p2[1] - p1[1]
        denom = n[0] * vx + n[1] * vy
        if abs(denom) < eps:
            continue
        t = (c - (n[0] * p1[0] + n[1] * p1[1])) / denom
        if t < -eps or t > 1 + eps:
            continue
        t = max(0.0, min(1.0, t))
        px = p1[0] + t * vx
        py = p1[1] + t * vy
        b_points.append((px, py))

    if len(b_points) < 2:
        return None

    # d 方向に並べて最小・最大
    s_vals = [d[0] * px + d[1] * py for (px, py) in b_points]
    s_min = min(s_vals)
    s_max = max(s_vals)
    if s_max - s_min < eps:
        return None

    # 対応する点
    idx_min = s_vals.index(s_min)
    idx_max = s_vals.index(s_max)
    p_min = b_points[idx_min]
    p_max = b_points[idx_max]
    L = s_max - s_min

    # 線分との交点（s の範囲内）
    hits = []
    for (x1, y1, x2, y2) in line_segments:
        p1 = (x1, y1)
        p2 = (x2, y2)
        vx = p2[0] - p1[0]
        vy = p2[1] - p1[1]
        denom = n[0] * vx + n[1] * vy
        if abs(denom) < eps:
            # 走査線と平行 (または同一直線) → 基本的に無視
            continue
        t = (c - (n[0] * p1[0] + n[1] * p1[1])) / denom
        if t < -eps or t > 1 + eps:
            continue
        t = max(0.0, min(1.0, t))
        px = p1[0] + t * vx
        py = p1[1] + t * vy
        s = d[0] * px + d[1] * py
        if s_min - eps <= s <= s_max + eps:
            hits.append((px, py))

    N = len(hits)
    density = N / L if L > 0 else float("nan")

    return {
        "pos": c,
        "length": L,
        "hits": N,
        "density": density,
        "boundary_points": [p_min, p_max],
        "hit_points": hits,
    }


# ---------------------- App ----------------------
class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Crackscanner V.1.0 ")
        self.geometry("1320x930")

        self.current_line_segments = []
        self.current_boundary_segments = []
        self.current_boundary_bbox = None

        # 単一走査線可視化用
        self.scan_y_val = None
        self.scan_hits = []
        self.scan_boundary_points = []
        self.scan_density = None
        self.scan_length = None

        # 水平複数走査線結果
        self.multi_scan_result = None  # dict or None
        self.multi_scan_draw = []      # 描画用: 各走査線の y, 交点など

        # 角度ごとの結果（0〜180°をステップ刻み）
        self.angle_results = {}    # angle_deg -> result dict (複数走査線結果)
        self.angle_draw = {}       # angle_deg -> draw list
        self.angle_conv = {}       # angle_deg -> (cum_lengths, cum_densities)
        self.angle_order = []      # タブ表示順
        self.current_angle = None  # 現在表示中の角度（deg）

        # Convergence plot用（別ウインドウ）
        self.conv_win = None
        self.conv_canvas = None
        self.conv_fig = None
        self.conv_ax = None

        # Rose diagram用（別ウインドウ）
        self.rose_win = None
        self.rose_canvas = None
        self.rose_fig = None
        self.rose_ax = None

        # ローズ図のExport CSV用データ（角度・密度など）
        self.rose_export_data = None

        # セグメント長さヒストグラム用
        self.segment_lengths = []
        self.length_win = None
        self.length_canvas = None
        self.length_fig = None
        self.length_ax = None

        # Controls (file)
        ctrl = tk.Frame(self)
        ctrl.pack(side=tk.TOP, fill=tk.X, padx=8, pady=6)
        tk.Button(ctrl, text="Open file...", command=self.open_file).pack(side=tk.LEFT)
        self.path_var = tk.StringVar(value="(no file selected)")
        tk.Label(ctrl, textvariable=self.path_var, anchor="w").pack(side=tk.LEFT, padx=10)

        # Controls (scale + axes)
        scale = tk.Frame(self)
        scale.pack(side=tk.TOP, fill=tk.X, padx=8, pady=4)
        # Ratio input
        self.units_val_var = tk.StringVar(value="1.0")
        self.mm_val_var = tk.StringVar(value="1.0")
        tk.Entry(scale, textvariable=self.units_val_var, width=7, justify="right").pack(side=tk.LEFT, padx=2)
        tk.Label(scale, text="Units =").pack(side=tk.LEFT)
        tk.Entry(scale, textvariable=self.mm_val_var, width=7, justify="right").pack(side=tk.LEFT, padx=2)
        tk.Label(scale, text="mm").pack(side=tk.LEFT, padx=(0, 12))

        # Projection choice
        tk.Label(scale, text="Axes:").pack(side=tk.LEFT, padx=(8, 0))
        self.proj_var = tk.StringVar(value="XY")
        tk.OptionMenu(scale, self.proj_var, "XY", "XZ", "YZ").pack(side=tk.LEFT, padx=2)

        # Scale bar length
        tk.Label(scale, text="Scale bar length").pack(side=tk.LEFT, padx=(12, 0))
        self.scale_mm_var = tk.StringVar(value="50")
        # ★ 0.5mm, 1.0mm を追加
        tk.OptionMenu(scale, self.scale_mm_var, "0.5", "1.0", "5", "10", "20", "50", "100", "200").pack(side=tk.LEFT, padx=2)
        tk.Label(scale, text="mm").pack(side=tk.LEFT)

        # Flip options: axis glyph only
        self.flip_x_var = tk.BooleanVar(value=False)
        self.flip_y_var = tk.BooleanVar(value=False)
        tk.Checkbutton(
            scale, text="Flip left-right (axes only)", variable=self.flip_x_var, command=self.redraw
        ).pack(side=tk.LEFT, padx=8)
        tk.Checkbutton(
            scale, text="Flip up-down (axes only)", variable=self.flip_y_var, command=self.redraw
        ).pack(side=tk.LEFT, padx=4)

        self.show_scale_var = tk.BooleanVar(value=True)
        tk.Checkbutton(
            scale, text="Show scale bar", variable=self.show_scale_var, command=self.redraw
        ).pack(side=tk.LEFT, padx=12)
        tk.Button(scale, text="Update", command=self.redraw).pack(side=tk.LEFT, padx=6)

        # Controls (single scanline)
        scan = tk.Frame(self)
        scan.pack(side=tk.TOP, fill=tk.X, padx=8, pady=4)
        tk.Label(scan, text="Single scanline y =").pack(side=tk.LEFT)
        self.scan_y_var = tk.StringVar(value="")
        tk.Entry(scan, textvariable=self.scan_y_var, width=10, justify="right").pack(side=tk.LEFT, padx=4)
        tk.Button(scan, text="Set center y of boundary", command=self.set_scanline_center).pack(side=tk.LEFT, padx=4)
        tk.Button(scan, text="Single scanline analysis", command=self.run_scanline_analysis).pack(side=tk.LEFT, padx=4)

        self.result_var = tk.StringVar(value="Single scanline analysis結果：-")
        tk.Label(scan, textvariable=self.result_var, anchor="w").pack(side=tk.LEFT, padx=10)

        # Controls (multi scanline)
        multi = tk.Frame(self)
        multi.pack(side=tk.TOP, fill=tk.X, padx=8, pady=4)

        # row 1: basic multi-scanline settings
        multi_row1 = tk.Frame(multi)
        multi_row1.pack(side=tk.TOP, fill=tk.X)

        tk.Label(multi_row1, text="Multi scanlines:").pack(side=tk.LEFT)

        self.multi_mode_var = tk.StringVar(value="even")
        tk.Radiobutton(
            multi_row1, text="Fixed number, evenly spaced", variable=self.multi_mode_var, value="even"
        ).pack(side=tk.LEFT, padx=4)
        tk.Radiobutton(
            multi_row1, text="Max number, sequential", variable=self.multi_mode_var, value="incremental"
        ).pack(side=tk.LEFT, padx=4)

        tk.Label(multi_row1, text="Number / max number =").pack(side=tk.LEFT, padx=(8, 2))
        self.multi_num_var = tk.StringVar(value="5")
        tk.Entry(multi_row1, textvariable=self.multi_num_var, width=7, justify="right").pack(side=tk.LEFT, padx=2)

        tk.Button(multi_row1, text="Multi scanline analysis (horizontal)", command=self.run_multi_scanlines).pack(side=tk.LEFT, padx=6)
        tk.Button(multi_row1, text="Export CSV", command=self.export_csv).pack(side=tk.LEFT, padx=4)

        # row 2: angle sweep settings
        multi_row2 = tk.Frame(multi)
        multi_row2.pack(side=tk.TOP, fill=tk.X, pady=(2, 0))

        tk.Label(multi_row2, text="Angle step").pack(side=tk.LEFT, padx=(8, 2))
        self.angle_step_var = tk.StringVar(value="15")
        tk.OptionMenu(multi_row2, self.angle_step_var, "15", "10").pack(side=tk.LEFT, padx=2)
        tk.Label(multi_row2, text="deg").pack(side=tk.LEFT)

        tk.Button(multi_row2, text="Angle sweep analysis", command=self.run_angle_sweep).pack(side=tk.LEFT, padx=6)

        self.multi_result_var = tk.StringVar(value="Multi scanlines: not computed")
        tk.Label(multi_row2, textvariable=self.multi_result_var, anchor="w").pack(side=tk.LEFT, padx=10)

# ローズ図用パネル
        rose_ctrl = tk.Frame(self)
        rose_ctrl.pack(side=tk.TOP, fill=tk.X, padx=8, pady=4)

        # row 1: L_total settings
        rose_row1 = tk.Frame(rose_ctrl)
        rose_row1.pack(side=tk.TOP, fill=tk.X)

        tk.Label(rose_row1, text="Rose: total scanline length L_total (units) =").pack(side=tk.LEFT)
        self.rose_L_var = tk.StringVar(value="")
        tk.Entry(rose_row1, textvariable=self.rose_L_var, width=10, justify="right").pack(side=tk.LEFT, padx=4)
        tk.Label(rose_row1, text="(Leave blank or 0 to use the full length for each angle)").pack(side=tk.LEFT, padx=4)

        # row 2: max radius and buttons
        rose_row2 = tk.Frame(rose_ctrl)
        rose_row2.pack(side=tk.TOP, fill=tk.X, pady=(2, 0))

        tk.Label(rose_row2, text="Max radius (count/mm, blank = auto) =").pack(side=tk.LEFT, padx=8)
        self.rose_rmax_var = tk.StringVar(value="")
        tk.Entry(rose_row2, textvariable=self.rose_rmax_var, width=10, justify="right").pack(side=tk.LEFT, padx=4)

        tk.Button(rose_row2, text="Show rose diagram", command=self.show_rose_diagram).pack(side=tk.LEFT, padx=8)
        tk.Button(rose_row2, text="Save rose SVG", command=self.save_rose_svg).pack(side=tk.LEFT, padx=4)

# セグメント長さヒストグラム用パネル
        length_ctrl = tk.Frame(self)
        length_ctrl.pack(side=tk.TOP, fill=tk.X, padx=8, pady=4)
        tk.Label(length_ctrl, text="Segment length:").pack(side=tk.LEFT)
        tk.Button(length_ctrl, text="Show length histogram", command=self.show_length_histogram).pack(side=tk.LEFT, padx=4)
        tk.Button(length_ctrl, text="Save histogram SVG", command=self.save_length_hist_svg).pack(side=tk.LEFT, padx=4)
        tk.Button(length_ctrl, text="lengthExport CSV", command=self.export_lengths_csv).pack(side=tk.LEFT, padx=4)
        
        # 角度タブ（メイン画面）
        angle_tab_frame = tk.Frame(self)
        angle_tab_frame.pack(side=tk.TOP, fill=tk.X, padx=8, pady=2)
        tk.Label(angle_tab_frame, text="Angle tabs:").pack(side=tk.LEFT)
        self.angle_notebook = ttk.Notebook(angle_tab_frame)
        self.angle_notebook.pack(side=tk.LEFT, fill=tk.X, expand=True)
        self.angle_notebook.bind("<<NotebookTabChanged>>", self.on_angle_tab_changed)

        # Figure for main drawing
        self.fig = Figure(figsize=(10.0, 7.4), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.ax.set_aspect("equal", adjustable="box")
        self.ax.axis("off")

        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        self.clear_plot()

    # ---------- Utility ----------
    def clear_plot(self):
        self.ax.clear()
        self.ax.set_aspect("equal", adjustable="box")
        self.ax.axis("off")
        self.canvas.draw_idle()

    def _mm_per_unit(self):
        try:
            u = float(self.units_val_var.get())
            m = float(self.mm_val_var.get())
            if u > 0:
                return m / u
        except Exception:
            pass
        return None

    # ---------- Convergence window ----------
    def open_convergence_window(self):
        # 既にウインドウがあれば前面に
        if self.conv_win is not None and self.conv_win.winfo_exists():
            self.conv_win.lift()
            return

        self.conv_win = tk.Toplevel(self)
        self.conv_win.title("Convergence plot")
        self.conv_win.geometry("600x400")

        self.conv_fig = Figure(figsize=(6, 4), dpi=100)
        self.conv_ax = self.conv_fig.add_subplot(111)
        self.conv_ax.set_title("Convergence plot")
        self.conv_ax.set_xlabel("Cumulative scanline length [units]")
        self.conv_ax.set_ylabel("Mean intersection density [count/unit]")
        self.conv_ax.grid(True, linestyle=":", linewidth=0.5)

        self.conv_canvas = FigureCanvasTkAgg(self.conv_fig, master=self.conv_win)
        self.conv_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.conv_canvas.draw()

    def update_convergence_plot(self, cum_lengths, cum_densities):
        # 別ウインドウを用意
        self.open_convergence_window()

        self.conv_ax.clear()
        self.conv_ax.plot(cum_lengths, cum_densities, marker="o", linestyle="-", linewidth=1.2)
        self.conv_ax.set_title("Convergence plot")
        self.conv_ax.set_xlabel("Cumulative scanline length [units]")
        self.conv_ax.set_ylabel("Mean intersection density [count/unit]")
        self.conv_ax.grid(True, linestyle=":", linewidth=0.5)
        self.conv_canvas.draw()

    # ---------- Rose window ----------
    def open_rose_window(self):
        if self.rose_win is not None and self.rose_win.winfo_exists():
            self.rose_win.lift()
            return

        self.rose_win = tk.Toplevel(self)
        self.rose_win.title("Rose diagram")
        self.rose_win.geometry("600x600")

        self.rose_fig = Figure(figsize=(6, 6), dpi=100)
        self.rose_ax = self.rose_fig.add_subplot(111, polar=True)

        self.rose_canvas = FigureCanvasTkAgg(self.rose_fig, master=self.rose_win)
        self.rose_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.rose_canvas.draw()


    def save_rose_svg(self):
        """
        表示中のRose diagramをSVG形式で保存する。
        先に「Show rose diagram」で図を描画してから実行することを想定。
        """
        if self.rose_fig is None:
            messagebox.showinfo("Information", 'The rose diagram has not been drawn yet. Please run "Show rose diagram" first.')
            return

        path = filedialog.asksaveasfilename(
            title="Save rose diagram as SVG",
            defaultextension=".svg",
            filetypes=[("SVG file", "*.svg"), ("All files", "*.*")],
        )
        if not path:
            return

        try:
            # SVG を保存
            self.rose_fig.savefig(path, format="svg")

            # ローズ図で使用したデータをCSVでも保存
            csv_msg = ""
            data = getattr(self, "rose_export_data", None)
            if data and isinstance(data, dict):
                try:
                    base, _ = os.path.splitext(path)
                    csv_path = base + "_rose_data.csv"
                    with open(csv_path, "w", newline="", encoding="utf-8") as f:
                        writer = csv.writer(f)
                        writer.writerow(["angle_deg", "density_per_mm", "L_total_used_unit", "mm_per_unit"])
                        angles = data.get("angles_deg") or []
                        dens = data.get("densities_per_mm") or []
                        L_target = data.get("L_target")
                        mm_per_unit = data.get("mm_per_unit")
                        for ang, d in zip(angles, dens):
                            writer.writerow([ang, d, L_target if L_target is not None else "", mm_per_unit])
                    csv_msg = f"\nRose-diagram data CSV: {csv_path}"
                except Exception as e_csv:
                    csv_msg = f"\n※ローズ図データのExport CSVに失敗しました: {e_csv}"

            messagebox.showinfo("Done", f"Save rose diagram as SVGしました:\n{path}{csv_msg}")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to save rose diagram:\n{e}")

    def show_rose_diagram(self):
        """
        Angle sweep analysis結果（self.angle_results）をもとに、
        各角度について「指定した総走査線長さ L_total」のところの平均密度を評価し、
        それを 0〜360°にミラー（対称）させたRose diagramを別ウインドウに表示。
        密度は本/単位 → 本/mm に変換して描画。
        """
        if not self.angle_results:
            messagebox.showinfo("Information", 'There are no angle-sweep multi-scanline results. Please run "Angle sweep analysis" first.')
            return

        mm_per_unit = self._mm_per_unit()
        if not mm_per_unit or mm_per_unit <= 0:
            messagebox.showerror("Error", "Scale (units→mm) is not set, so the rose diagram cannot be drawn in count/mm. Please set the conversion between units and mm.")
            return

        # L_total の取得（単位長さベース）
        L_target = None
        txt = self.rose_L_var.get().strip()
        if txt:
            try:
                val = float(txt)
                if val > 0:
                    L_target = val
            except Exception:
                messagebox.showerror("Error", "Please enter a numeric value for L_total.")
                return

        # 各角度で密度を決定（本/mm）
        base_angles = sorted(self.angle_results.keys())  # 0〜<180
        base_values = []

        for ang in base_angles:
            r = self.angle_results[ang]
            cumL = r["cum_lengths"]       # [単位]
            cumD = r["cum_densities"]     # [本/単位]
            if not cumL:
                continue

            if L_target is None:
                # その角度で利用可能な全長さ → 最終点の平均密度
                dens_unit = cumD[-1]
            else:
                # cumL が L_target を超える最初の点を採用
                dens_unit = cumD[-1]
                for L_val, d_val in zip(cumL, cumD):
                    if L_val >= L_target:
                        dens_unit = d_val
                        break

            dens_mm = dens_unit / mm_per_unit  # 本/mm
            base_values.append(dens_mm)

        if not base_values:
            messagebox.showinfo("Information", "No valid rose-diagram data were obtained.")
            return

        # CSV 出力用にデータを保存（角度[deg], 密度[本/mm]）
        self.rose_export_data = {
            "mm_per_unit": mm_per_unit,
            "L_target": L_target,
            "angles_deg": base_angles,
            "densities_per_mm": base_values,
        }

        # 角度間隔（バーの幅）を推定
        if len(base_angles) > 1:
            diffs = [b - a for a, b in zip(base_angles[:-1], base_angles[1:])]
            step_deg = min(diffs)
        else:
            step_deg = 10.0
        width_rad = math.radians(step_deg)

        # シンメトリーを考慮して 0〜360°にミラー
        thetas = []
        vals = []
        for ang, val in zip(base_angles, base_values):
            th = math.radians(ang)
            thetas.append(th)
            vals.append(val)
            th2 = math.radians(ang + 180.0)
            thetas.append(th2)
            vals.append(val)

        # ローズウインドウを用意して描画
        self.open_rose_window()

        self.rose_ax.clear()
        # 0°を上（北）に、時計回り
        self.rose_ax.set_theta_zero_location("N")
        self.rose_ax.set_theta_direction(-1)

        self.rose_ax.bar(
            thetas,
            vals,
            width=width_rad,
            bottom=0.0,
            align="center",
            edgecolor="black",
            linewidth=0.8,
            alpha=0.8,
        )

        if L_target is None:
            title = "Rose diagram（/mm）"
        else:
            title = f"Rose diagram（L_total ≈ {L_target:.2f} 単位までの平均密度: 本/mm）"
        self.rose_ax.set_title(title)

        # --- ローズ図上に座標軸ラベル（x, y など）を追加 ---
        try:
            proj = (self.proj_var.get() or "XY").upper()
        except Exception:
            proj = "XY"
        if proj == "XY":
            hx_label = "x"
            vy_label = "y"
        elif proj == "XZ":
            hx_label = "x"
            vy_label = "z"
        else:
            hx_label = "y"
            vy_label = "z"

        # 左右・上下反転設定に応じて方位（度）を決定
        # 0°/180° に水平方向（x または y）、90°/270° に鉛直方向（y または z）を割り当てる。
        ang_x = 0.0 if not self.flip_x_var.get() else 180.0
        ang_y = 90.0 if not self.flip_y_var.get() else 270.0

        th_x = math.radians(ang_x)
        th_y = math.radians(ang_y)

        # 半径方向の位置（少し外側に）※ vals はバーの高さ配列
        r_max_data = max(vals) if vals else 1.0

        # 最大半径（軸スケール）のオプション指定
        rmax_spec = None
        try:
            if hasattr(self, "rose_rmax_var"):
                txt_r = self.rose_rmax_var.get().strip()
                if txt_r:
                    v = float(txt_r)
                    if v <= 0:
                        raise ValueError
                    rmax_spec = v
        except ValueError:
            messagebox.showerror("Error", "Please enter a positive numeric value for the rose-diagram maximum radius.")
            return
        except Exception:
            # 入力が取得できない場合は自動スケールに任せる
            rmax_spec = None

        if rmax_spec is not None:
            r_max = rmax_spec
            # 軸スケールを指定値に固定（デフォルトは自動スケールのまま）
            self.rose_ax.set_ylim(0.0, r_max)
        else:
            r_max = r_max_data

        r_label = r_max * 1.05

        # ラベル描画
        self.rose_ax.text(
            th_x,
            r_label,
            hx_label,
            ha="center",
            va="center",
            fontsize=12,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=1.0),
        )
        self.rose_ax.text(
            th_y,
            r_label,
            vy_label,
            ha="center",
            va="center",
            fontsize=12,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=1.0),
        )

        self.rose_canvas.draw()

    # ---------- セグメント長さヒストグラムウインドウ ----------
    def open_length_window(self):
        if self.length_win is not None and self.length_win.winfo_exists():
            self.length_win.lift()
            return

        self.length_win = tk.Toplevel(self)
        self.length_win.title("Segment-length histogram")
        self.length_win.geometry("600x400")

        self.length_fig = Figure(figsize=(6, 4), dpi=100)
        self.length_ax = self.length_fig.add_subplot(111)

        self.length_canvas = FigureCanvasTkAgg(self.length_fig, master=self.length_win)
        self.length_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.length_canvas.draw()


    def save_length_hist_svg(self):
        """
        表示中のSegment-length histogramをSVG形式で保存する。
        先に「Show length histogram」で図を描画してから実行することを想定。
        """
        if self.length_fig is None:
            messagebox.showinfo("Information", 'The length histogram has not been drawn yet. Please run "Show length histogram" first.')
            return

        path = filedialog.asksaveasfilename(
            title="Save length histogram as SVG",
            defaultextension=".svg",
            filetypes=[("SVG file", "*.svg"), ("All files", "*.*")],
        )
        if not path:
            return

        try:
            self.length_fig.savefig(path, format="svg")
            messagebox.showinfo("Done", f"Save length histogram as SVGしました:\n{path}")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to save length histogram:\n{e}")

    def show_length_histogram(self):
        # 線分の端点同士の距離（長さ）のヒストグラムを別ウインドウに表示。
        # 縦軸は規格化した頻度（確率）で、全ビンの高さの総和が 1 になるようにする。
        if not self.segment_lengths:
            messagebox.showinfo("Information", "No segment-length data. Please load a file first.")
            return

        mm_per_unit = self._mm_per_unit()
        if not mm_per_unit or mm_per_unit <= 0:
            messagebox.showerror("Error", "Scale (units→mm) is not set, so the length histogram in mm cannot be drawn. Please set the conversion between units and mm.")
            return

        self.open_length_window()

        self.length_ax.clear()

        # NumPy のヒストグラム機能でビンを自動決定し、規格化する（mm単位）
        lengths_arr = np.array(self.segment_lengths, dtype=float)
        lengths_mm = lengths_arr * mm_per_unit
        counts, bin_edges = np.histogram(lengths_mm, bins="auto")
        total = counts.sum()
        if total <= 0:
            messagebox.showinfo("Information", "No valid data for histogram.")
            return

        probs = counts / total  # 各ビンの確率（合計 1）
        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
        bin_widths = bin_edges[1:] - bin_edges[:-1]

        self.length_ax.bar(
            bin_centers,
            probs,
            width=bin_widths,
            edgecolor="black",
            linewidth=0.8,
        )
        self.length_ax.set_title("Segment-length histogram（unit:mm）")
        self.length_ax.set_xlabel("Segment length [mm]")
        self.length_ax.set_ylabel("Normalized frequency (probability)")
        self.length_ax.grid(True, linestyle=":", linewidth=0.5)

        self.length_canvas.draw()

    def export_lengths_csv(self):
        # 線分長さデータを CSV 出力する（単位＆mm 両方）。
        if not self.segment_lengths:
            messagebox.showinfo("Information", "No segment-length data. Please load a file first.")
            return

        mm_per_unit = self._mm_per_unit()
        if not mm_per_unit or mm_per_unit <= 0:
            messagebox.showerror("Error", "Scale (units→mm) is not set, so CSV in mm cannot be exported. Please set the conversion between units and mm.")
            return

        path = filedialog.asksaveasfilename(
            title="Select destination for segment-length CSV",
            defaultextension=".csv",
            filetypes=[("CSV file", "*.csv"), ("All files", "*.*")],
        )
        if not path:
            return

        try:
            with open(path, "w", newline="", encoding="utf-8") as f:
                writer = csv.writer(f)
                # index, 単位長さ, mm長さ
                writer.writerow(["index", "length_units", "length_mm"])
                for i, L_units in enumerate(self.segment_lengths, start=1):
                    L_mm = L_units * mm_per_unit
                    writer.writerow([i, f"{L_units:.6f}", f"{L_mm:.6f}"])
            messagebox.showinfo("Done", f"Exported segment-length CSV:\n{path}")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to export segment-length CSV:\n{e}")

    # ---------- Drawing ----------
    def draw_all(self, line_segments, boundary_segments, boundary_bbox):
        self.ax.clear()
        self.ax.set_aspect("equal", adjustable="box")
        self.ax.axis("off")

        xs, ys = [], []
        for (x1, y1, x2, y2) in line_segments:
            self.ax.plot([x1, x2], [y1, y2], "-", linewidth=0.8, color="black")
            xs.extend([x1, x2])
            ys.extend([y1, y2])

        for (x1, y1, x2, y2) in boundary_segments:
            self.ax.plot([x1, x2], [y1, y2], "-", linewidth=1.2, color="red")
            xs.extend([x1, x2])
            ys.extend([y1, y2])

        if not xs:
            self.canvas.draw_idle()
            return

        xmin, xmax = min(xs), max(xs)
        ymin, ymax = min(ys), max(ys)
        dx = xmax - xmin if xmax > xmin else 1.0
        dy = ymax - ymin if ymax > ymin else 1.0
        mx = dx * 0.06
        my = dy * 0.06
        base_xlim = (xmin - mx, xmax + mx)
        base_ylim = (ymin - my, ymax + my)

        # Boundary bbox fallback
        if boundary_bbox:
            xmin_b, xmax_b, ymin_b, ymax_b = boundary_bbox
        else:
            xmin_b, xmax_b, ymin_b, ymax_b = xmin, xmax, ymin, ymax

        # --- Scale bar (bottom-right OUTSIDE) ---
        bar_geom = None
        mm_per_unit = self._mm_per_unit()
        if self.show_scale_var.get() and mm_per_unit and mm_per_unit > 0:
            try:
                bar_mm = float(self.scale_mm_var.get())
                if bar_mm > 0:
                    bar_len_units = bar_mm / mm_per_unit
                    gap_y = max(0.05 * dy, 1.0)
                    gap_x = max(0.04 * dx, 1.0)
                    x0 = xmax_b + gap_x
                    y0 = ymax_b + gap_y
                    x1 = x0 + bar_len_units
                    y1 = y0
                    cap = max(0.012 * dy, 0.5)
                    voff = max(0.035 * dy, 1.0)
                    label = f"{float(bar_mm):.1f} mm" if abs(bar_mm - round(bar_mm)) > 1e-6 else f"{int(round(bar_mm))} mm"
                    bar_geom = (x0, y0, x1, y1, cap, voff, label)

                    pad_x = 0.06 * dx
                    pad_y_top = 0.06 * dy
                    ex_xmin = min(base_xlim[0], x0 - cap * 0.5 - pad_x)
                    ex_xmax = max(base_xlim[1], x1 + cap * 0.5 + pad_x)
                    ex_ymin = min(base_ylim[0], y0 - cap - voff - pad_y_top)
                    ex_ymax = max(base_ylim[1], y0 + cap + pad_y_top)
                    base_xlim = (ex_xmin, ex_xmax)
                    base_ylim = (ex_ymin, ex_ymax)
            except Exception:
                bar_geom = None

        # --- Mini axis glyph (bottom-left OUTSIDE, 原点固定) ---
        axis_geom = None
        try:
            L = max(0.12 * dy, 2.0)
            gap_y2 = max(0.06 * dy, 1.0)
            gap_x2 = max(0.06 * dx, 1.0)
            xg = xmin_b - gap_x2
            yg = ymax_b + gap_y2
            axis_geom = (xg, yg, L)

            pad_x2 = 0.08 * dx + 0.20 * L
            pad_y2 = 0.08 * dy + 0.20 * L
            ex_xmin = min(base_xlim[0], xg - pad_x2)
            ex_ymax = max(base_ylim[1], yg + pad_y2)
            base_xlim = (ex_xmin, base_xlim[1])
            base_ylim = (base_ylim[0], ex_ymax)
        except Exception:
            axis_geom = None

        # 図の向きは固定（flip_x/flip_yでは xlim/ylim をいじらない）
        self.ax.set_xlim(*base_xlim)
        self.ax.set_ylim(*base_ylim)
        self.ax.invert_yaxis()

        # スケールバー描画
        if bar_geom:
            x0, y0, x1, y1, cap, voff, label = bar_geom
            self.ax.plot([x0, x1], [y0, y1], "-", linewidth=2.0, color="black")
            self.ax.plot([x0, x0], [y0 - cap, y0 + cap], "-", linewidth=2.0, color="black")
            self.ax.plot([x1, x1], [y1 - cap, y1 + cap], "-", linewidth=2.0, color="black")
            self.ax.text(
                (x0 + x1) / 2.0,
                y0 - voff,
                label,
                ha="center",
                va="top",
                fontsize=11,
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=1.5),
            )

        # ミニ座標軸
        if axis_geom:
            xg, yg, L = axis_geom
            hx_dir = -1 if self.flip_x_var.get() else 1
            vy_dir = 1 if self.flip_y_var.get() else -1

            self.ax.plot([xg, xg + L * hx_dir], [yg, yg], "-", linewidth=2.0, color="black")
            self.ax.plot([xg, xg], [yg, yg + L * vy_dir], "-", linewidth=2.0, color="black")

            ah = max(0.025 * L, 0.8)

            xt = xg + L * hx_dir
            xb = xt - ah * hx_dir
            self.ax.plot([xt, xb], [yg, yg - ah], "-", linewidth=2.0, color="black")
            self.ax.plot([xt, xb], [yg, yg + ah], "-", linewidth=2.0, color="black")

            yt = yg + L * vy_dir
            yb = yt - ah * vy_dir
            self.ax.plot([xg, xg - ah], [yt, yb], "-", linewidth=2.0, color="black")
            self.ax.plot([xg, xg + ah], [yt, yb], "-", linewidth=2.0, color="black")

            proj = (self.proj_var.get() or "XY").upper()
            if proj == "XY":
                hx, vy = "x", "y"
            elif proj == "XZ":
                hx, vy = "x", "z"
            else:
                hx, vy = "y", "z"

            hx_label_x = xt + 0.10 * L * hx_dir
            hx_label_y = yg
            vy_label_x = xg
            vy_label_y = yt + 0.12 * L * vy_dir

            self.ax.text(
                hx_label_x,
                hx_label_y,
                hx,
                ha="left" if hx_dir > 0 else "right",
                va="center",
                fontsize=12,
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=1.5),
            )
            self.ax.text(
                vy_label_x,
                vy_label_y,
                vy,
                ha="center",
                va="top" if vy_dir < 0 else "bottom",
                fontsize=12,
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=1.5),
            )

        # --- 単一走査線の描画（水平のみ） ---
        if self.scan_y_val is not None and self.current_boundary_bbox is not None:
            y0 = self.scan_y_val
            self.ax.plot(
                [base_xlim[0], base_xlim[1]],
                [y0, y0],
                "--",
                linewidth=1.5,
                color="blue",
            )
            for (xb, yb) in self.scan_boundary_points:
                self.ax.plot(xb, yb, "o", color="blue", markersize=5)
            for (xh, yh) in self.scan_hits:
                self.ax.plot(xh, yh, "x", color="blue", markersize=5, mew=1.5)

        # --- 複数走査線の描画（水平 / 角度付き） ---
        if self.multi_scan_draw:
            for sd in self.multi_scan_draw:
                # 水平の場合は 'y' を、角度付きの場合は 'line_points' を持つようにする
                if "y" in sd:
                    y0 = sd["y"]
                    # 図全体を横切る破線
                    self.ax.plot(
                        [base_xlim[0], base_xlim[1]],
                        [y0, y0],
                        "--",
                        linewidth=1.0,
                        color="blue",
                        alpha=0.5,
                    )
                elif "line_points" in sd:
                    (px0, py0), (px1, py1) = sd["line_points"]
                    self.ax.plot(
                        [px0, px1],
                        [py0, py1],
                        "--",
                        linewidth=1.0,
                        color="blue",
                        alpha=0.5,
                    )

                # 境界との交点
                for (xb, yb) in sd["boundary_points"]:
                    self.ax.plot(xb, yb, "o", color="blue", markersize=3, alpha=0.7)
                # 線分との交点
                for (xh, yh) in sd["hit_points"]:
                    self.ax.plot(xh, yh, "x", color="blue", markersize=3, mew=1.0, alpha=0.7)

        self.canvas.draw_idle()

    # ---------- Scanline helpers (horizontal) ----------
    def set_scanline_center(self):
        if not self.current_boundary_bbox:
            messagebox.showinfo("Information", "No boundary found. Please load a file first.")
            return
        xmin_b, xmax_b, ymin_b, ymax_b = self.current_boundary_bbox
        yc = 0.5 * (ymin_b + ymax_b)
        self.scan_y_var.set(f"{yc:.3f}")

    def analyze_scanline(self, y0):
        """
        従来の水平Single scanline analysis。
        """
        if not self.current_boundary_bbox:
            return None

        # 走査線と境界の交点
        b_intersections = []
        for (x1, y1, x2, y2) in self.current_boundary_segments:
            x = _intersect_segment_with_horizontal(x1, y1, x2, y2, y0)
            if x is not None:
                b_intersections.append((x, y0))
        b_intersections.sort(key=lambda p: p[0])

        if len(b_intersections) < 2:
            return None

        left = b_intersections[0]
        right = b_intersections[-1]
        L = abs(right[0] - left[0])
        if L <= 0:
            return None

        x_min, x_max = min(left[0], right[0]), max(left[0], right[0])
        hits = []
        for (x1, y1, x2, y2) in self.current_line_segments:
            x = _intersect_segment_with_horizontal(x1, y1, x2, y2, y0)
            if x is not None and x_min <= x <= x_max:
                hits.append((x, y0))

        N = len(hits)
        density = N / L

        return {
            "y": y0,
            "length": L,
            "hits": N,
            "density": density,
            "boundary_points": [left, right],
            "hit_points": hits,
        }

    def run_scanline_analysis(self):
        if not self.current_boundary_bbox:
            messagebox.showerror("Error", "Cannot run scanline analysis because the boundary (red box) was not found.")
            return
        try:
            y0 = float(self.scan_y_var.get())
        except Exception:
            messagebox.showerror("Error", "Please enter a numeric value for scanline y.")
            return

        res = self.analyze_scanline(y0)
        if res is None:
            self.scan_y_val = y0
            self.scan_boundary_points = []
            self.scan_hits = []
            self.scan_length = None
            self.scan_density = None
            self.result_var.set(f"Single scanline analysis結果：y={y0:.3f} で有効な区間が得られませんでした")
            self.redraw()
            return

        self.scan_y_val = res["y"]
        self.scan_boundary_points = res["boundary_points"]
        self.scan_hits = res["hit_points"]
        self.scan_length = res["length"]
        self.scan_density = res["density"]

        mm_per_unit = self._mm_per_unit()
        if mm_per_unit and mm_per_unit > 0:
            L_mm = res["length"] * mm_per_unit
            rho_mm = res["hits"] / L_mm if L_mm > 0 else float("nan")
            rho_cm = rho_mm * 10.0
            txt = (
                f"y={res['y']:.3f}: N={res['hits']}, L={res['length']:.3f} units, "
                f"ρ={res['density']:.4f} count/unit "
                f"(L={L_mm:.2f} mm, ρ={rho_mm:.4f} count/mm, {rho_cm:.4f} count/cm)"
            )
        else:
            txt = (
                f"y={res['y']:.3f}: N={res['hits']}, L={res['length']:.3f} units, "
                f"ρ={res['density']:.4f} count/unit"
            )

        self.result_var.set("Single scanline analysis結果：" + txt)
        self.redraw()

    # ---------- Multi scanlines (horizontal, existing) ----------
    def run_multi_scanlines(self):
        """
        従来の「Multi scanline analysis (horizontal)」。
        """
        if not self.current_boundary_bbox:
            messagebox.showerror("Error", "Cannot run multi-scanline analysis because the boundary (red box) was not found.")
            return
        try:
            n = int(self.multi_num_var.get())
            if n <= 0:
                raise ValueError
        except Exception:
            messagebox.showerror("Error", "Please enter an integer ≥ 1 for number / max number.")
            return

        mode = self.multi_mode_var.get()  # "even" or "incremental"
        xmin_b, xmax_b, ymin_b, ymax_b = self.current_boundary_bbox
        H = ymax_b - ymin_b
        if H <= 0:
            messagebox.showerror("Error", "Boundary height is zero.")
            return

        ys = []
        for i in range(1, n + 1):
            t = i / (n + 1)
            y0 = ymin_b + t * H
            ys.append(y0)

        lengths = []
        hits_list = []
        densities = []
        valid_indices = []
        multi_draw = []

        for idx, y0 in enumerate(ys, start=1):
            res = self.analyze_scanline(y0)
            if res is None:
                continue
            lengths.append(res["length"])
            hits_list.append(res["hits"])
            densities.append(res["density"])
            valid_indices.append(idx)
            multi_draw.append(
                {
                    "y": res["y"],
                    "boundary_points": res["boundary_points"],
                    "hit_points": res["hit_points"],
                }
            )

        if not lengths:
            self.multi_scan_result = None
            self.multi_scan_draw = []
            self.multi_result_var.set("Multi scanlines: no valid scanlines were obtained")

            if self.conv_ax is not None:
                self.conv_ax.clear()
                self.conv_ax.set_title("Convergence plot")
                self.conv_ax.set_xlabel("Cumulative scanline length [units]")
                self.conv_ax.set_ylabel("Mean intersection density [count/unit]")
                self.conv_ax.grid(True, linestyle=":", linewidth=0.5)
                if self.conv_canvas is not None:
                    self.conv_canvas.draw()
            self.canvas.draw_idle()
            return

        cum_lengths = []
        cum_hits = []
        cum_densities = []
        sL = 0.0
        sN = 0
        for L, N in zip(lengths, hits_list):
            sL += L
            sN += N
            cum_lengths.append(sL)
            cum_hits.append(sN)
            cum_densities.append(sN / sL if sL > 0 else float("nan"))

        self.multi_scan_result = {
            "mode": mode,
            "direction": "horizontal_y",
            "indices": valid_indices,
            "ys": [ys[i - 1] for i in valid_indices],
            "lengths": lengths,
            "hits": hits_list,
            "densities": densities,
            "cum_lengths": cum_lengths,
            "cum_hits": cum_hits,
            "cum_densities": cum_densities,
        }
        self.multi_scan_draw = multi_draw

        total_L = cum_lengths[-1]
        total_N = cum_hits[-1]
        mean_rho = cum_densities[-1]
        self.multi_result_var.set(
            f"Multi scanlines (horizontal): valid count {len(lengths)}, total length {total_L:.3f} units, "
            f"総交点数 {total_N}, Mean density {mean_rho:.4f} count/unit"
        )

        self.update_convergence_plot(cum_lengths, cum_densities)
        self.redraw()

    # ---------- Angle sweep (0〜180°を step°刻み) ----------
    def run_angle_sweep(self):
        if not self.current_boundary_bbox:
            messagebox.showerror("Error", "Cannot run angle-sweep multi-scanline analysis because the boundary (red box) was not found.")
            return
        try:
            n = int(self.multi_num_var.get())
            if n <= 0:
                raise ValueError
        except Exception:
            messagebox.showerror("Error", "Please enter an integer ≥ 1 for number / max number.")
            return

        try:
            step = int(self.angle_step_var.get())
            if step <= 0 or step > 180:
                step = 15
        except Exception:
            step = 15

        mode = self.multi_mode_var.get()  # "even" or "incremental"
        # 角度リスト: 0,step,...,<180
        angles = list(range(0, 180, step))

        if not self.current_boundary_segments:
            messagebox.showerror("Error", "No boundary segments.")
            return

        # 境界の端点を列挙（n・p の min/max 用）
        boundary_points = []
        for (x1, y1, x2, y2) in self.current_boundary_segments:
            boundary_points.append((x1, y1))
            boundary_points.append((x2, y2))

        angle_results = {}
        angle_draw = {}
        angle_conv = {}
        angle_order = []

        for ang in angles:
            theta = math.radians(ang)
            d = (math.cos(theta), math.sin(theta))
            nvec = (-math.sin(theta), math.cos(theta))

            # nvec・p の min/max
            vals = [nvec[0] * px + nvec[1] * py for (px, py) in boundary_points]
            if not vals:
                continue
            cmin = min(vals)
            cmax = max(vals)
            if cmax - cmin <= 1e-9:
                continue

            # 等間隔に nvec・p = c を取る（端を避ける）
            cs = []
            for i in range(1, n + 1):
                t = i / (n + 1)
                c = cmin + t * (cmax - cmin)
                cs.append(c)

            lengths = []
            hits_list = []
            densities = []
            valid_indices = []
            multi_draw_ang = []

            for idx, c in enumerate(cs, start=1):
                res = analyze_scanline_general(
                    self.current_line_segments,
                    self.current_boundary_segments,
                    d,
                    nvec,
                    c,
                )
                if res is None:
                    continue
                lengths.append(res["length"])
                hits_list.append(res["hits"])
                densities.append(res["density"])
                valid_indices.append(idx)

                # 線分としての描画用に、boundary_points の 2 点を使用
                p_min, p_max = res["boundary_points"]
                multi_draw_ang.append(
                    {
                        "line_points": [p_min, p_max],
                        "boundary_points": res["boundary_points"],
                        "hit_points": res["hit_points"],
                    }
                )

            if not lengths:
                continue

            cum_lengths = []
            cum_hits = []
            cum_densities = []
            sL = 0.0
            sN = 0
            for L, N in zip(lengths, hits_list):
                sL += L
                sN += N
                cum_lengths.append(sL)
                cum_hits.append(sN)
                cum_densities.append(sN / sL if sL > 0 else float("nan"))

            angle_results[ang] = {
                "mode": mode,
                "direction": f"angle_{ang}deg",
                "indices": valid_indices,
                "ys": cs,  # 一般角度では n・p=c を「位置」として格納
                "lengths": lengths,
                "hits": hits_list,
                "densities": densities,
                "cum_lengths": cum_lengths,
                "cum_hits": cum_hits,
                "cum_densities": cum_densities,
            }
            angle_draw[ang] = multi_draw_ang
            angle_conv[ang] = (cum_lengths, cum_densities)
            angle_order.append(ang)

        if not angle_order:
            messagebox.showinfo("Information", "No valid scanlines were obtained for any angle.")
            return

        # 保存
        self.angle_results = angle_results
        self.angle_draw = angle_draw
        self.angle_conv = angle_conv
        self.angle_order = sorted(angle_order)
        self.current_angle = self.angle_order[0]

        # 角度タブを作り直し
        for tab_id in self.angle_notebook.tabs():
            self.angle_notebook.forget(tab_id)
        for ang in self.angle_order:
            frame = ttk.Frame(self.angle_notebook)
            self.angle_notebook.add(frame, text=f"{ang}°")

        # 最初の角度を表示
        ang0 = self.current_angle
        self.multi_scan_draw = self.angle_draw.get(ang0, [])
        self.multi_scan_result = self.angle_results.get(ang0)

        r = self.multi_scan_result
        if r:
            total_L = r["cum_lengths"][-1]
            total_N = r["cum_hits"][-1]
            mean_rho = r["cum_densities"][-1]
            self.multi_result_var.set(
                f"Angle {ang0}°: valid count {len(r['lengths'])}, "
                f"Total length {total_L:.3f} units, total intersections {total_N}, "
                f"Mean density {mean_rho:.4f} count/unit"
            )

            self.update_convergence_plot(r["cum_lengths"], r["cum_densities"])

        self.redraw()

    def on_angle_tab_changed(self, event):
        # 角度タブ切り替え時 → multi_scan_draw / Convergence plotを切り替え
        if not self.angle_order:
            return
        cur_idx = self.angle_notebook.index(self.angle_notebook.select())
        if cur_idx < 0 or cur_idx >= len(self.angle_order):
            return
        ang = self.angle_order[cur_idx]
        self.current_angle = ang

        if ang in self.angle_draw:
            self.multi_scan_draw = self.angle_draw[ang]
        if ang in self.angle_results:
            self.multi_scan_result = self.angle_results[ang]
            r = self.multi_scan_result
            total_L = r["cum_lengths"][-1]
            total_N = r["cum_hits"][-1]
            mean_rho = r["cum_densities"][-1]
            self.multi_result_var.set(
                f"Angle {ang}°: valid count {len(r['lengths'])}, "
                f"Total length {total_L:.3f} units, total intersections {total_N}, "
                f"Mean density {mean_rho:.4f} count/unit"
            )
            self.update_convergence_plot(r["cum_lengths"], r["cum_densities"])

        self.redraw()

    # ---------- CSV export ----------
    def export_csv(self):
        """
        CSV 出力:
        - 角度付き解析結果（self.angle_results）があれば、すべての角度・走査線を 1 ファイルに出力
        - なければ、従来の self.multi_scan_result（水平方向のみ）を出力
        出力は単位と mm の両方を含む。
        """
        if self.angle_results:
            mode = "angle"
        elif self.multi_scan_result:
            mode = "horizontal"
        else:
            messagebox.showinfo("Information", "There are no multi-scanline analysis results to export. Please run the analysis first.")
            return

        mm_per_unit = self._mm_per_unit()
        if not mm_per_unit or mm_per_unit <= 0:
            messagebox.showerror("Error", "Scale (units→mm) is not set, so CSV in mm cannot be exported. Please set the conversion between units and mm.")
            return

        path = filedialog.asksaveasfilename(
            title="CSV fileとして保存",
            defaultextension=".csv",
            filetypes=[("CSV file", "*.csv"), ("All files", "*.*")],
        )
        if not path:
            return

        try:
            with open(path, "w", newline="", encoding="utf-8") as f:
                writer = csv.writer(f)
                # 共通ヘッダ
                writer.writerow(
                    [
                        "angle_deg",              # 角度（水平のみの場合は 0）
                        "mode",
                        "direction",
                        "index",
                        "y_or_pos_units",        # 水平: y[単位], 角度付き: n・p=c [単位]
                        "length_units",
                        "length_mm",
                        "hits",
                        "density_per_unit",
                        "density_per_mm",
                        "cum_length_units",
                        "cum_length_mm",
                        "cum_hits",
                        "cum_density_per_unit",
                        "cum_density_per_mm",
                    ]
                )

                if mode == "angle":
                    for ang in sorted(self.angle_results.keys()):
                        r = self.angle_results[ang]
                        for idx, pos, L, N, d, cL, cN, cD in zip(
                            r["indices"],
                            r["ys"],
                            r["lengths"],
                            r["hits"],
                            r["densities"],
                            r["cum_lengths"],
                            r["cum_hits"],
                            r["cum_densities"],
                        ):
                            L_mm = L * mm_per_unit
                            d_mm = d / mm_per_unit if mm_per_unit > 0 else float("nan")
                            cL_mm = cL * mm_per_unit
                            cD_mm = cD / mm_per_unit if mm_per_unit > 0 else float("nan")
                            writer.writerow(
                                [
                                    ang,
                                    r["mode"],
                                    r["direction"],
                                    idx,
                                    f"{pos:.6f}",
                                    f"{L:.6f}",
                                    f"{L_mm:.6f}",
                                    N,
                                    f"{d:.8f}",
                                    f"{d_mm:.8f}",
                                    f"{cL:.6f}",
                                    f"{cL_mm:.6f}",
                                    cN,
                                    f"{cD:.8f}",
                                    f"{cD_mm:.8f}",
                                ]
                            )
                else:
                    r = self.multi_scan_result
                    for idx, y0, L, N, d, cL, cN, cD in zip(
                        r["indices"],
                        r["ys"],
                        r["lengths"],
                        r["hits"],
                        r["densities"],
                        r["cum_lengths"],
                        r["cum_hits"],
                        r["cum_densities"],
                    ):
                        L_mm = L * mm_per_unit
                        d_mm = d / mm_per_unit if mm_per_unit > 0 else float("nan")
                        cL_mm = cL * mm_per_unit
                        cD_mm = cD / mm_per_unit if mm_per_unit > 0 else float("nan")
                        writer.writerow(
                            [
                                0,  # 水平は angle=0 として出力
                                r["mode"],
                                r["direction"],
                                idx,
                                f"{y0:.6f}",
                                f"{L:.6f}",
                                f"{L_mm:.6f}",
                                N,
                                f"{d:.8f}",
                                f"{d_mm:.8f}",
                                f"{cL:.6f}",
                                f"{cL_mm:.6f}",
                                cN,
                                f"{cD:.8f}",
                                f"{cD_mm:.8f}",
                            ]
                        )

            messagebox.showinfo("Done", f"Saved CSV:\n{path}")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to export CSV:\n{e}")

    # ---------- IO ----------
    def open_file(self):
        path = filedialog.askopenfilename(
            title="SVG fileを選択",
            filetypes=[("SVG / Text", "*.svg *.txt"), ("All files", "*.*")],
        )
        if not path:
            return
        self.path_var.set(path)
        try:
            with open(path, "r", encoding="utf-8", errors="ignore") as f:
                text = f.read()
            ls, bs, bb = extract_all(text)
        except Exception as e:
            messagebox.showerror("解析Error", f"Failed to load/analyze file:\n{e}")
            return

        self.current_line_segments = ls
        self.current_boundary_segments = bs
        self.current_boundary_bbox = bb

        # クラック長さ配列（ヒストグラム用）をUpdate
        # ここでの「クラック」は SVG 中の line / polyline 要素を 1 本として扱い，
        # polyline の長さは最初と最後の点の距離とする。
        self.segment_lengths = compute_crack_lengths(text)

        # 単一・複数走査線Informationをリセット
        self.scan_y_val = None
        self.scan_boundary_points = []
        self.scan_hits = []
        self.scan_length = None
        self.scan_density = None
        self.result_var.set("Single scanline analysis結果：-")

        self.multi_scan_result = None
        self.multi_scan_draw = []
        self.multi_result_var.set("Multi scanlines: not computed")

        # 角度付きもリセット
        self.angle_results = {}
        self.angle_draw = {}
        self.angle_conv = {}
        self.angle_order = {}
        self.current_angle = None
        for tab_id in self.angle_notebook.tabs():
            self.angle_notebook.forget(tab_id)

        # 収束ウインドウが開いていればクリア
        if self.conv_ax is not None:
            self.conv_ax.clear()
            self.conv_ax.set_title("Convergence plot")
            self.conv_ax.set_xlabel("Cumulative scanline length [units]")
            self.conv_ax.set_ylabel("Mean intersection density [count/unit]")
            self.conv_ax.grid(True, linestyle=":", linewidth=0.5)
            if self.conv_canvas is not None:
                self.conv_canvas.draw()

        self.draw_all(ls, bs, bb)

    def redraw(self):
        self.draw_all(self.current_line_segments, self.current_boundary_segments, self.current_boundary_bbox)


def main():
    app = App()
    if len(sys.argv) > 1 and os.path.isfile(sys.argv[1]):
        try:
            with open(sys.argv[1], "r", encoding="utf-8", errors="ignore") as f:
                text = f.read()
            ls, bs, bb = extract_all(text)
            app.path_var.set(sys.argv[1])
            app.current_line_segments = ls
            app.current_boundary_segments = bs
            app.current_boundary_bbox = bb

            # クラック長さ配列（ヒストグラム用）をUpdate
            # ここでの「クラック」は SVG 中の line / polyline 要素を 1 本として扱い，
            # polyline の長さは最初と最後の点の距離とする。
            app.segment_lengths = compute_crack_lengths(text)

            app.draw_all(ls, bs, bb)
        except Exception as e:
            messagebox.showerror("解析Error", f"Failed to analyze {sys.argv[1]}:\n{e}")
    app.mainloop()


if __name__ == "__main__":
    main()
