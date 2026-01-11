from PIL import Image, ImageDraw, ImageFont
import os

# -----------------------------------------------------------
# CONFIGURATION
# -----------------------------------------------------------
WIDTH = 400
HEIGHT = 6000

PIXELS_PER_DEGREE = 20
DEG_MIN = -90
DEG_MAX = 90

# E-bar geometry
LONG_SEG = 120
SHORT_SEG = 60
GAP = 20

# HUD color
COLOR = (0, 255, 0, 255)

# --- Updated font and thickness ---
FONT_SIZE = 40
FONT_PATH = r"C:\Windows\Fonts\arial.ttf"  # bold font
TEXT_STROKE = 2

MAJOR_WIDTH = 6
MINOR_WIDTH = 4
ZERO_WIDTH = 6
# -----------------------------------------------------------

# Create transparent image
img = Image.new("RGBA", (WIDTH, HEIGHT), (0, 0, 0, 0))
draw = ImageDraw.Draw(img)
font = ImageFont.truetype(FONT_PATH, FONT_SIZE)

center_y = HEIGHT // 2

for angle in range(DEG_MIN, DEG_MAX + 1, 5):
    y = center_y - angle * PIXELS_PER_DEGREE

    if y < 0 or y >= HEIGHT:
        continue

    # Special zero-degree tick
    if angle == 0:
        draw.line((100, y, 300, y), fill=COLOR, width=ZERO_WIDTH)
        continue

    # Major tick at every 10 degrees
    if angle % 10 == 0:
        # Left long segment
        draw.line((40, y, 40 + LONG_SEG, y), fill=COLOR, width=MAJOR_WIDTH)
        # Center short segment
        x_center = 40 + LONG_SEG + GAP
        draw.line((x_center, y, x_center + SHORT_SEG, y), fill=COLOR, width=MAJOR_WIDTH)
        # Right long segment
        x_right = x_center + SHORT_SEG + GAP
        draw.line((x_right, y, x_right + LONG_SEG, y), fill=COLOR, width=MAJOR_WIDTH)
    else:
        # Minor tick every 5 degrees
        draw.line((140, y, 260, y), fill=COLOR, width=MINOR_WIDTH)

    # Draw the angle number for every 5 degrees
    num = str(abs(angle))  # Using absolute like F-16 style
    bbox = draw.textbbox((0, 0), num, font=font, stroke_width=TEXT_STROKE)
    text_w = bbox[2] - bbox[0]
    text_h = bbox[3] - bbox[1]

    # Position text above for positive, below for negative
    if angle > 0:
        txt_y = y - 50
    else:
        txt_y = y + 20

    draw.text(
        (WIDTH // 2 - text_w // 2, txt_y),
        num,
        fill=COLOR,
        font=font,
        stroke_fill=COLOR,
        stroke_width=TEXT_STROKE
    )

# Save output
script_dir = os.path.dirname(os.path.abspath(__file__))
output_path = os.path.join(script_dir, "pitch_bars.png")
img.save(output_path)
print("Saved:", output_path)
