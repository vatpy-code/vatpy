import os
import tkinter as tk
from tkinter import ttk
from PIL import Image, ImageTk


class ImageViewerApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Vatpy Image Viewer")

        # Label for dropdown
        self.label = ttk.Label(root, text="Select an image:")
        self.label.pack(pady=5)

        # StringVar for dropdown
        self.image_var = tk.StringVar()

        # Dropdown menu for image selection
        self.image_dropdown = ttk.Combobox(root, textvariable=self.image_var)
        self.image_dropdown.pack(pady=5)

        # Button to refresh the list of images
        self.refresh_button = ttk.Button(root, text="Refresh List",
                                         command=self.update_image_list)
        self.refresh_button.pack(pady=5)

        # Button to display the selected image
        self.display_button = ttk.Button(root, text="Display Image",
                                         command=self.display_image)
        self.display_button.pack(pady=10)

        # Label to display the image
        self.image_label = ttk.Label(root)
        self.image_label.pack(pady=10)

        # Initialize the list of images
        self.update_image_list()

    def update_image_list(self):
        # Update the list of available image files in the current directory:
        self.image_files = [
            f for f in os.listdir()
            if f.lower().endswith(('.png', '.jpg', '.jpeg', '.gif'))
        ]
        self.image_dropdown['values'] = self.image_files
        if self.image_files:
            self.image_var.set(self.image_files[0])

    def display_image(self):
        selected_image = self.image_var.get()
        if selected_image in self.image_files:
            try:
                img = Image.open(selected_image)
                img.thumbnail((700, 700))
                img_tk = ImageTk.PhotoImage(img)
                self.image_label.config(image=img_tk)
                self.image_label.image = img_tk
            except Exception as e:
                self.image_label.config(text=f"Error: {e}")
        else:
            self.image_label.config(text="Please select a valid image.")


if __name__ == "__main__":
    root = tk.Tk()
    app = ImageViewerApp(root)
    root.mainloop()
