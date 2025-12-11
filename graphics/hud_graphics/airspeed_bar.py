from PIL import Image, ImageDraw, ImageFont
import os

# -----------------------------------------------------------
# CONFIGURATION
# -----------------------------------------------------------
WIDTH = 400
HEIGHT = 24000

PIXELS_PER_UNIT = 10
MIN_VAL = 0
MAX_VAL = 1600

MINOR_STEP = 10
MAJOR_STEP = 50

LONG_SEG = 200
SHORT_SEG = 120
COLOR = (0, 255, 0, 255)

FONT_SIZE = 40   # bigger text
font = ImageFont.truetype(r"C:\Windows\Fonts\arial.ttf", FONT_SIZE)  # bold font
TEXT_STROKE = 2  # outline thickness for text

MINOR_WIDTH = 4  # thicker minor ticks
MAJOR_WIDTH = 6  # thicker major ticks
# -----------------------------------------------------------

img = Image.new("RGBA", (WIDTH, HEIGHT), (0, 0, 0, 0))
draw = ImageDraw.Draw(img)

def val_to_y(v):
    return HEIGHT - int(v * PIXELS_PER_UNIT)

for v in range(MIN_VAL, MAX_VAL + MINOR_STEP, MINOR_STEP):

    y = val_to_y(v)

    if not (0 <= y < HEIGHT):
        continue

    # -------- Major tick every 50 --------
    if v % MAJOR_STEP == 0:

        draw.line((40, y, 40 + LONG_SEG, y), fill=COLOR, width=MAJOR_WIDTH)
        draw.line((WIDTH - 40 - LONG_SEG, y, WIDTH - 40, y), fill=COLOR, width=MAJOR_WIDTH)

        # Label
        num = str(v)
        bbox = draw.textbbox((0, 0), num, font=font, stroke_width=TEXT_STROKE)
        text_w = bbox[2] - bbox[0]
        text_h = bbox[3] - bbox[1]

        # Zero label special placement
        if v == 0:
            draw.text(
                (WIDTH // 2 - text_w // 2, y + 6),
                num, fill=COLOR, font=font, stroke_fill=COLOR, stroke_width=TEXT_STROKE
            )
        else:
            draw.text(
                (WIDTH // 2 - text_w // 2, y - text_h - 10),
                num, fill=COLOR, font=font, stroke_fill=COLOR, stroke_width=TEXT_STROKE
            )

    # -------- Minor ticks every 10 --------
    else:
        draw.line(
            (120, y, WIDTH - 120, y),
            fill=COLOR,
            width=MINOR_WIDTH
        )

script_dir = os.path.dirname(os.path.abspath(__file__))
output_path = os.path.join(script_dir, "airspeed_bars.png")
img.save(output_path)

print("Saved:", output_path)
