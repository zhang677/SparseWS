def pngs_to_pdf(folder: str):
  """
  Input: Folder name with images
  Output: Combine all images into one pdf. Each 12 images are on one page.
  """
  import os
  from PIL import Image
  from fpdf import FPDF

  # get all png files in folder
  files = [f for f in os.listdir(folder) if f.endswith('.png')]
  files.sort()

  # create pdf
  pdf = FPDF()
  pdf.set_auto_page_break(0)

  # add all images to pdf
  for i, file in enumerate(files):
    if i % 12 == 0:
      pdf.add_page()
    pdf.image(os.path.join(folder, file), x=(i % 3) * 65, y=int(i / 3) % 4 * 65, w=65, h=65)

  # save pdf
  pdf.output(os.path.join(folder, 'output.pdf'))

if __name__ == '__main__':
  pngs_to_pdf('/home/nfs_data/zhanggh/SparseWS/data/monitor')