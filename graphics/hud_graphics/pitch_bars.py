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

FONT_SIZE = 30
# FONT_PATH = "Arial.ttf"  # Change to your preferred HUD font
# -----------------------------------------------------------

# Create transparent image
img = Image.new("RGBA", (WIDTH, HEIGHT), (0, 0, 0, 0))
draw = ImageDraw.Draw(img)
font = ImageFont.truetype(r"C:\Windows\Fonts\arial.ttf", 32)

center_y = HEIGHT // 2

for angle in range(DEG_MIN, DEG_MAX + 1, 5):

    y = center_y - angle * PIXELS_PER_DEGREE

    if y < 0 or y >= HEIGHT:
        continue

    # Major tick at every 10 degrees (E-bar)
    if angle % 10 == 0 and angle != 0:
        # Positive pitch: E-bar points upward
        if angle > 0:
            bar_y = y - 10
            txt_y = y - 50
        # Negative pitch: E-bar points downward
        else:
            bar_y = y + 10
            txt_y = y + 20

        # Left long segment
        draw.line((40, bar_y, 40 + LONG_SEG, bar_y), fill=COLOR, width=3)
        # Center short segment
        x_center = 40 + LONG_SEG + GAP
        draw.line((x_center, bar_y, x_center + SHORT_SEG, bar_y), fill=COLOR, width=3)
        # Right long segment
        x_right = x_center + SHORT_SEG + GAP
        draw.line((x_right, bar_y, x_right + LONG_SEG, bar_y), fill=COLOR, width=3)

        # Draw angle number
        num = str(abs(angle))  # F-16 uses absolute for negatives
        bbox = draw.textbbox((0, 0), num, font=font)
        text_w = bbox[2] - bbox[0]
        text_h = bbox[3] - bbox[1]        
        draw.text((WIDTH // 2 - text_w // 2, txt_y), num, fill=COLOR, font=font)

    # Minor tick every 5 degrees (simple line)
    elif angle != 0:
        draw.line(
            (140, y, 260, y),
            fill=COLOR,
            width=2
        )

# Save output
script_dir = os.path.dirname(os.path.abspath(__file__))
output_path = os.path.join(script_dir, "f16_pitch_ladder.png")

img.save(output_path)
print("Saved:", output_path)