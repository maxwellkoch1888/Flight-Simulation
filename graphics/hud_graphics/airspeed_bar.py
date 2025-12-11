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
MAJOR_STEP = 100

LONG_SEG = 200
SHORT_SEG = 120
COLOR = (0, 255, 0, 255)

FONT_SIZE = 30
font = ImageFont.truetype(r"C:\Windows\Fonts\arial.ttf", FONT_SIZE)
# -----------------------------------------------------------

img = Image.new("RGBA", (WIDTH, HEIGHT), (0, 0, 0, 0))
draw = ImageDraw.Draw(img)

def val_to_y(v):
    return HEIGHT - int(v * PIXELS_PER_UNIT)

for v in range(MIN_VAL, MAX_VAL + MINOR_STEP, MINOR_STEP):

    y = val_to_y(v)

    if not (0 <= y < HEIGHT):
        continue

    # -------- Major tick every 100 --------
    if v % MAJOR_STEP == 0:

        draw.line((40, y, 40 + LONG_SEG, y), fill=COLOR, width=3)
        draw.line((WIDTH - 40 - LONG_SEG, y, WIDTH - 40, y), fill=COLOR, width=3)

        # Label
        num = str(v)
        bbox = draw.textbbox((0, 0), num, font=font)
        text_w = bbox[2] - bbox[0]
        text_h = bbox[3] - bbox[1]

        # Special placement for 0 label
        if v == 0:
            # put label *below* the line
            draw.text((WIDTH // 2 - text_w // 2, y + 4),
                      num, fill=COLOR, font=font)
        else:
            # normal placement (above line)
            draw.text((WIDTH // 2 - text_w // 2, y - text_h - 8),
                      num, fill=COLOR, font=font)

    else:
        draw.line(
            (120, y, WIDTH - 120, y),
            fill=COLOR,
            width=2
        )

script_dir = os.path.dirname(os.path.abspath(__file__))
output_path = os.path.join(script_dir, "airspeed_bars.png")
img.save(output_path)

print("Saved:", output_path)
