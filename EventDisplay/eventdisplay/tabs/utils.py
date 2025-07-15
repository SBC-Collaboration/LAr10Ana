from tkinter import Label
from PIL import Image, ImageTk
from pathlib import Path

def show_not_found_image(parent, message):
    """
    Display an error message and fallback image inside a Tkinter container.
    """
    image_path = Path(__file__).parent.parent.parent / "notfound.jpeg"
    image_size = (300, 200)
    font_settings = ("TkDefaultFont", 14)

    # Show error message
    Label(parent, text=message, fg='red', font=font_settings).grid(row=0, column=0, sticky='NW', padx=10, pady=10)

    # Load and display image if it exists
    if image_path.exists():
        try:
            image = Image.open(image_path).resize(image_size)
            photo = ImageTk.PhotoImage(image)

            label = Label(parent, image=photo)
            label.image = photo  # Prevent garbage collection
            label.grid(row=1, column=0, padx=10, pady=10)
        except Exception as e:
            Label(parent, text=f"[Image error: {e}]", fg='gray', font=font_settings).grid(row=1, column=0, sticky='NW')
    else:
        Label(parent, text="notfound.jpeg is missing", fg='gray', font=font_settings).grid(row=1, column=0, sticky='NW')
