from PIL import Image, ImageDraw, ImageFont
import os

# -----------------------------------------------------------
# CONFIGURATION
# -----------------------------------------------------------
WIDTH = 400
HEIGHT = 24000

PIXELS_PER_FOOT = 0.5
ALT_MIN = 0
ALT_MAX = 100000

BOTTOM_MARGIN_FT = 50
MAPPED_MIN = ALT_MIN - BOTTOM_MARGIN_FT

MINOR_STEP = 100
MAJOR_STEP = 500

LONG_SEG = 200
SHORT_SEG = 120
COLOR = (0, 255, 0, 255)

# --- Updated font and thickness settings ---
FONT_SIZE = 40
font = ImageFont.truetype(r"C:\Windows\Fonts\arial.ttf", FONT_SIZE)

TEXT_STROKE = 2

MAJOR_WIDTH = 6   # thicker major ticks
MINOR_WIDTH = 4   # thicker minor ticks
# -----------------------------------------------------------

# Create transparent image
img = Image.new("RGBA", (WIDTH, HEIGHT), (0, 0, 0, 0))
draw = ImageDraw.Draw(img)

def alt_to_y(alt_ft):
    return HEIGHT - int((alt_ft - MAPPED_MIN) * PIXELS_PER_FOOT)

for alt in range(ALT_MIN, ALT_MAX + MINOR_STEP, MINOR_STEP):

    y = alt_to_y(alt)

    if y < 0 or y >= HEIGHT:
        continue

    # -------- Major tick every 500 ft --------
    if alt % MAJOR_STEP == 0:

        draw.line((40, y, 40 + LONG_SEG, y), fill=COLOR, width=MAJOR_WIDTH)
        draw.line((WIDTH - 40 - LONG_SEG, y, WIDTH - 40, y), fill=COLOR, width=MAJOR_WIDTH)

        num = str(alt)
        bbox = draw.textbbox((0, 0), num, font=font, stroke_width=TEXT_STROKE)
        text_w = bbox[2] - bbox[0]
        text_h = bbox[3] - bbox[1]

        draw.text(
            (WIDTH // 2 - text_w // 2, y - text_h - 10),
            num,
            fill=COLOR,
            font=font,
            stroke_fill=COLOR,
            stroke_width=TEXT_STROKE
        )

    # -------- Minor tick every 100 ft --------
    else:
        draw.line(
            (120, y, WIDTH - 120, y),
            fill=COLOR,
            width=MINOR_WIDTH
        )

script_dir = os.path.dirname(os.path.abspath(__file__))
output_path = os.path.join(script_dir, "altitude_bars.png")
img.save(output_path)

print("Saved:", output_path)
