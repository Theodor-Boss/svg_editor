# lige en lille bemærkning

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk


class Point:
    def __init__(self, x, y, tag=None):
        self.x = x
        self.y = y
        self.tag = tag  # Optional tag (e.g., color or label)


class LineSegment:
    def __init__(self, start, end, tag=None):
        self.x0, self.y0 = start
        self.x1, self.y1 = end
        self.tag = tag  # Optional tag (e.g., color or label)

    def line_segment(self, t):
        x = self.x0 + t * (self.x1 - self.x0)
        y = self.y0 + t * (self.y1 - self.y0)
        return x, y


class Cubic_bezier:
    def __init__(self, start, a, b, end, tag=None):
        self.x0, self.y0 = start
        self.x1, self.y1 = a
        self.x2, self.y2 = b
        self.x3, self.y3 = end
        self.tag = tag

    def cubic_bezier(self, t):
        x = (1-t)**3 * self.x0 + 3*(1-t)**2 * t * self.x1 + 3*(1-t) * t**2 * self.x2 + t**3 * self.x3
        y = (1-t)**3 * self.y0 + 3*(1-t)**2 * t * self.y1 + 3*(1-t) * t**2 * self.y2 + t**3 * self.y3
        return x, y


class Curves:
    def __init__(self):
        self.curves = []

    def add_curve(self, curve):
        self.curves.append(curve)

    def __iter__(self):
        return iter(self.curves)


def extract_elements(instructions, tag):
    subunits = []
    subunit = ''
    for i, character in enumerate(instructions):
        if character in ['C', 'L', 'M', 'Q'] and len(subunit) > 1:
            subunits.append(subunit)
            subunit = ''
        subunit += character
    start_x, start_y = list(map(float, subunits[0][1:].split()))
    kurver = Curves()
    # kurver.add_curve(Point(start_x, -start_y, tag=tag))
    for thing in subunits[1:]:
        numbers = list(map(float, thing[1:].split()))
        if thing[0] == 'C':
            new_spline = Cubic_bezier((start_x, -start_y), (numbers[0], -numbers[1]), (numbers[2], -numbers[3]), (numbers[4], -numbers[5]), tag=tag)
            kurver.add_curve(new_spline)

            start_x = numbers[-2]
            start_y = numbers[-1]

        elif thing[0] == 'L':
            new_line = LineSegment((start_x, -start_y), (numbers[0], -numbers[1]), tag=tag)
            kurver.add_curve(new_line)
            start_x = numbers[-2]
            start_y = numbers[-1]

        elif thing[0] == 'Q':
            start_x = numbers[-2]
            start_y = numbers[-1]
    return kurver


def line_is_relevant(key_word, line):
    for i in range(len(line) - len(key_word) + 1):
        piece = line[i:i + len(key_word)]
        if piece == key_word:
            return True
    return False


def trim_line(key_word1, key_word2, line):
    i = 0
    piece = ""
    while not piece == key_word1:
        i += 1
        piece = line[i:i + len(key_word1)]
    i += len(key_word1) - 1

    j = len(line) - len(key_word2) - 1
    piece = ""
    while not piece == key_word2:
        j -= 1
        piece = line[j:j + len(key_word2)]
    return line[i:j]


class SVGEditorApp:
    def __init__(self, root):
        self.root = root
        self.root.title("SVG Curve Selector")

        # Properly handle window close event
        self.root.protocol("WM_DELETE_WINDOW", self.on_closing)

        # Store curves and selections
        self.all_curves = Curves()
        self.selected_curves = set()  # This set will store indices of deleted curves

        # Mouse-related attributes
        self.start_x, self.start_y = None, None
        self.rect_id = None

        # Create and set up the plot figure and axes
        self.fig, self.ax = plt.subplots()
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.root)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        # Bind mouse wheel event for zooming
        self.canvas.mpl_connect("scroll_event", self.zoom)

        # Bind mouse events
        # self.canvas.get_tk_widget().bind("<Button-1>", self.start_delete_on_drag)
        # self.canvas.get_tk_widget().bind("<B1-Motion>", self.delete_on_drag)
        self.canvas.get_tk_widget().bind("<Button-3>", self.start_select_rectangle)
        self.canvas.get_tk_widget().bind("<B3-Motion>", self.draw_rectangle)
        self.canvas.get_tk_widget().bind("<ButtonRelease-3>", self.delete_in_rectangle)

        # Create Listbox to show curve elements
        self.curve_listbox = tk.Listbox(self.root, selectmode=tk.MULTIPLE)
        # self.curve_listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=1)

        # Create buttons for deletion and updating
        self.delete_button = tk.Button(self.root, text="Delete Selected Curves", command=self.delete_selected_curves)
        # self.delete_button.pack(side=tk.BOTTOM, fill=tk.X)

        self.initial_xlim = None
        self.initial_ylim = None
        self.y_max = None  # Store the maximum y value for flipping coordinates

        # Load the SVG file and populate the curves
        self.load_svg("fish_image.svg")
        # self.load_svg("dont_give_up_on_your.svg")

    def on_closing(self):
        """Handle cleanup when the window is closed."""
        plt.close(self.fig)  # Close the Matplotlib figure
        self.root.destroy()  # Destroy the Tkinter window and exit the application

    def load_svg(self, file_path):
        with open(file_path, "r") as fish:
            lines = fish.readlines()

        colors = ['green', 'blue', 'red', 'orange', 'purple', 'yellow']
        current_color = 0
        kriterie = r'd="M'
        kriterie2 = r'Z"'
        for line in lines:
            if line_is_relevant(kriterie, line):
                current_color = (current_color + 1) % len(colors)
                instructions = trim_line(kriterie, kriterie2, line)
                my_curves = extract_elements(instructions, tag=colors[current_color])
                for curve in my_curves:
                    self.all_curves.add_curve(curve)

        # Find the maximum y value
        self.y_max = max(max(curve.y0, curve.y1, curve.y2, curve.y3) if isinstance(curve, Cubic_bezier) else
                         max(curve.y0, curve.y1) if isinstance(curve, LineSegment) else
                         curve.y
                         for curve in self.all_curves)

        # Populate the listbox
        for idx, curve in enumerate(self.all_curves):
            self.curve_listbox.insert(tk.END, f"Curve {idx+1} ({type(curve).__name__})")

        # Plot all curves initially
        self.update_plot(True)

    def update_plot(self, reset_view=False):
        current_xlim = self.ax.get_xlim()
        current_ylim = self.ax.get_ylim()

        self.ax.clear()

        # Iterate over all curves to plot them
        for idx, element in enumerate(self.all_curves):
            if idx in self.selected_curves:
                continue  # Skip any deleted curves

            # Plot the curve based on its type
            if isinstance(element, Cubic_bezier):
                t_vals = np.linspace(0, 1, 100)
                x_curve, y_curve = element.cubic_bezier(t_vals)
                self.ax.plot(x_curve, y_curve, color=element.tag)
            elif isinstance(element, LineSegment):
                t_vals = np.linspace(0, 1, 100)
                x_curve, y_curve = element.line_segment(t_vals)
                self.ax.plot(x_curve, y_curve, color=element.tag)
            elif isinstance(element, Point):
                self.ax.scatter(element.x, element.y, c=element.tag)

        if reset_view or self.initial_xlim is None:
            self.ax.set_aspect('equal', adjustable='datalim')
            self.initial_xlim = self.ax.get_xlim()
            self.initial_ylim = self.ax.get_ylim()
        else:
            self.ax.set_xlim(current_xlim)
            self.ax.set_ylim(current_ylim)

        self.canvas.draw()

    def delete_selected_curves(self):
        # Get selected indices from the listbox
        selected_indices = list(self.curve_listbox.curselection())

        # Add selected indices to the set of deleted curves
        self.selected_curves.update(selected_indices)

        # Remove selected curves from the list
        for idx in sorted(selected_indices, reverse=True):
            del self.all_curves.curves[idx]
            self.curve_listbox.delete(idx)

        # Update plot after deletion
        self.update_plot(False)

    def start_delete_on_drag(self, event):
        """Start the deletion of curves as the user drags the left mouse button."""
        self.delete_on_drag(event)

    def delete_on_drag(self, event):
        """Delete any curve touched by the left mouse drag."""
        x, y = event.x, event.y
        x_data, y_data = self.ax.transData.inverted().transform((x, y))

        for idx, element in enumerate(self.all_curves):
            if idx in self.selected_curves:
                continue
            if isinstance(element, Point):
                if np.hypot(element.x - x_data, element.y - y_data) < 0.05:
                    self.selected_curves.add(idx)
            elif isinstance(element, LineSegment):
                x0, y0, x1, y1 = element.x0, element.y0, element.x1, element.y1
                if min(x0, x1) <= x_data <= max(x0, x1) and min(y0, y1) <= y_data <= max(y0, y1):
                    self.selected_curves.add(idx)
            elif isinstance(element, Cubic_bezier):
                t_vals = np.linspace(0, 1, 100)
                x_curve, y_curve = element.cubic_bezier(t_vals)
                if any(np.hypot(x_curve - x_data, y_curve - y_data) < 0.05):
                    self.selected_curves.add(idx)

        self.update_plot(False)

    def start_select_rectangle(self, event):
        """Start drawing the selection rectangle for deletion."""
        self.start_x, self.start_y = event.x, event.y

    def draw_rectangle(self, event):
        """Draw the selection rectangle as the user drags the right mouse button."""
        if self.rect_id:
            self.canvas.get_tk_widget().delete(self.rect_id)
        x1, y1 = self.start_x, self.start_y
        x2, y2 = event.x, event.y
        self.rect_id = self.canvas.get_tk_widget().create_rectangle(x1, y1, x2, y2, outline="red", dash=(2, 2))

    def delete_in_rectangle(self, event):
        """Delete any curves touched by the selection rectangle."""
        if self.rect_id:
            self.canvas.get_tk_widget().delete(self.rect_id)
            self.rect_id = None

        x1, y1 = self.ax.transData.inverted().transform((self.start_x, self.canvas.get_tk_widget().winfo_height() - self.start_y))
        x2, y2 = self.ax.transData.inverted().transform((event.x, self.canvas.get_tk_widget().winfo_height() - event.y))

        x_min, x_max = min(x1, x2), max(x1, x2)
        y_min, y_max = min(y1, y2), max(y1, y2)

        for idx, element in enumerate(self.all_curves):
            if idx in self.selected_curves:
                continue
            if isinstance(element, Point):
                if x_min <= element.x <= x_max and y_min <= element.y <= y_max:
                    self.selected_curves.add(idx)
            elif isinstance(element, LineSegment):
                x0, y0, x1, y1 = element.x0, element.y0, element.x1, element.y1
                if (x_min <= x0 <= x_max and y_min <= y0 <= y_max) or (x_min <= x1 <= x_max and y_min <= y1 <= y_max):
                    self.selected_curves.add(idx)
            elif isinstance(element, Cubic_bezier):
                t_vals = np.linspace(0, 1, 100)
                x_curve, y_curve = element.cubic_bezier(t_vals)
                if any((x_min <= x <= x_max and y_min <= y <= y_max) for x, y in zip(x_curve, y_curve)):
                    self.selected_curves.add(idx)

        self.update_plot(False)

    def zoom(self, event):
        """Zoom in or out centered on the mouse position."""
        base_scale = 1.1  # Zoom factor

        # Get the current axis limits
        cur_xlim = self.ax.get_xlim()
        cur_ylim = self.ax.get_ylim()

        xdata = event.xdata  # Get mouse x position in data coordinates
        ydata = event.ydata  # Get mouse y position in data coordinates

        # Calculate zooming direction
        if event.button == 'up':  # Scroll up to zoom in
            scale_factor = 1 / base_scale
        elif event.button == 'down':  # Scroll down to zoom out
            scale_factor = base_scale
        else:
            scale_factor = 1

        # Calculate new limits
        new_width = (cur_xlim[1] - cur_xlim[0]) * scale_factor
        new_height = (cur_ylim[1] - cur_ylim[0]) * scale_factor

        # Calculate new axis limits, keeping the mouse cursor in the same position
        relx = (cur_xlim[1] - xdata) / (cur_xlim[1] - cur_xlim[0])
        rely = (cur_ylim[1] - ydata) / (cur_ylim[1] - cur_ylim[0])

        self.ax.set_xlim([xdata - new_width * (1 - relx), xdata + new_width * relx])
        self.ax.set_ylim([ydata - new_height * (1 - rely), ydata + new_height * rely])

        self.canvas.draw_idle()  # Redraw the canvas


if __name__ == "__main__":
    root = tk.Tk()
    app = SVGEditorApp(root)
    root.mainloop()
