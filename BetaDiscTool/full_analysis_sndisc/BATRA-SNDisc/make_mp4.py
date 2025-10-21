import sys
import glob
from PIL import Image

def create_gif(image_pattern, output_path, duration=400):
  image_files = sorted(glob.glob(image_pattern))
  if not image_files:
    print(f"No images found matching {image_pattern}")
    sys.exit(1)

  images = [Image.open(img) for img in image_files]

  images[0].save(output_path, save_all=True, append_images=images[1:], duration=duration, loop=0)

  print(f"Saved GIF to {output_path}")

if __name__ == "__main__":
  if len(sys.argv) != 3:
    print("Usage: python make_gif.py '<image_pattern>' <output_path.gif>")
    sys.exit(1)

  image_pattern = sys.argv[1]
  output_path = sys.argv[2]
  create_gif(image_pattern, output_path)
