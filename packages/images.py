from PIL import Image
from PIL import ImageDraw, ImageFont
import glob

"""input"""
input_path = "/Users/kotaro/PycharmProjects/iden/inputs/遺伝画像/"
output_path = "/Users/kotaro/PycharmProjects/iden/outputs/images/"
font_path = "/Users/kotaro/PycharmProjects/iden/inputs/NotoSansJP-Black.ttf"
samples = ["LacZ", "CG4319", "CG12284"]
green = {i: {} for i in samples}
red = {i: {} for i in samples}

sample_stages = {"LacZ": [5, 8, 9, 11, 17, 13],
                 "CG4319": [8, 9, 13],
                 "CG12284": [8, 11, 5],
                 }

"""definition"""
font = ImageFont.truetype(font_path, 30)
font2 = ImageFont.truetype(font_path, 40)


def get_concat_3(lst, title):
    name1, name2, name3 = lst[0][0], lst[1][0], lst[2][0]
    im1, im2, im3 = lst[0][1], lst[1][1], lst[2][1]
    dst = Image.new("RGB", (im1.width + im2.width + im3.width, im1.height))
    dst.paste(im1, (0, 0))
    dst.paste(im2, (im1.width, 0))
    dst.paste(im3, (im1.width + im2.width, 0))
    draw = ImageDraw.Draw(dst)
    draw.text((0, im1.height * 0), title, 'white', font=font2)
    draw.text((0, im1.height * 0.9), name1, 'white', font=font)
    draw.text((im1.width, im1.height * 0.9), name2, 'white', font=font)
    draw.text((im1.width + im2.width, im1.height * 0.9), name3, 'white', font=font)
    return dst


def get_concat_6(lst, title):
    name1, name2, name3, name4, name5, name6 = lst[0][0], lst[1][0], lst[2][0], lst[3][0], lst[4][0], lst[5][0]
    im1, im2, im3, im4, im5, im6 = lst[0][1], lst[1][1], lst[2][1], lst[3][1], lst[4][1], lst[5][1]
    dst = Image.new("RGB", (im1.width * 3, im1.height * 2))
    dst.paste(im1, (0, 0))
    dst.paste(im2, (im1.width, 0))
    dst.paste(im3, (im1.width + im2.width, 0))
    dst.paste(im4, (0, im1.height))
    dst.paste(im5, (im1.width, im1.height))
    dst.paste(im6, (im1.width + im2.width, im1.height))
    draw = ImageDraw.Draw(dst)
    draw.text((0, im1.height * 0), title, 'white', font=font2)
    draw.text((0, im1.height * 0.9), name1, 'white', font=font)
    draw.text((im1.width, im1.height * 0.9), name2, 'white', font=font)
    draw.text((im1.width + im2.width, im1.height * 0.9), name3, 'white', font=font)
    draw.text((0, im1.height * 1.9), name4, 'white', font=font)
    draw.text((im1.width, im1.height * 1.9), name5, 'white', font=font)
    draw.text((im1.width + im2.width, im1.height * 1.9), name6, 'white', font=font)
    return dst

def images():
    """get_data"""
    stages = {}
    for key in sample_stages:
        stages[key] = ["STAGE" + str(i) for i in sample_stages[key]]

    order = {}
    for key, value in sample_stages.items():
        order[key] = ["STAGE" + str(i) for i in sorted(value)]

    for i in glob.glob(input_path + "ImageJ_緑/*/*.jpg"):
        green[i.split("/")[8]][stages[i.split("/")[8]][int(i[-15]) - 1]] = Image.open(i)

    for i in glob.glob(input_path + "ImageJ_赤/*/*.jpg"):
        red[i.split("/")[8]][stages[i.split("/")[8]][int(i[-17]) - 1]] = Image.open(i)

    """same_stages"""
    same_stages_img = {}
    same_stages = set(stages["LacZ"]) & set(stages["CG4319"]) & set(stages["CG12284"])

    for key in same_stages:
        same_stages_img[key] = []
        for i in green:
            same_stages_img[key].append([i, green[i][key]])

    same_stages_img_red = {}
    for key in same_stages:
        same_stages_img_red[key] = []
        for i in green:
            same_stages_img_red[key].append([i, red[i][key]])

    for key in same_stages_img:
        get_concat_3(same_stages_img[key], key).save(output_path + "merged " + str(key) + ".jpg")

    for key in same_stages_img_red:
        get_concat_3(same_stages_img_red[key], key).save(output_path + "merged " + str(key) + "(red).jpg")

    """same_samples"""
    same_samples_img = {"LacZ": [],
                        "CG4319": [],
                        "CG12284": []
                        }

    for sample in order:
        for i in order[sample]:
            same_samples_img[sample].append([i, green[sample][i]])

    for key in same_samples_img:
        if len(same_samples_img[key]) == 6:
            get_concat_6(same_samples_img[key], key).save(output_path + "merged " + str(key) + ".jpg")
        if len(same_samples_img[key]) == 3:
            get_concat_3(same_samples_img[key], key).save(output_path + "merged " + str(key) + ".jpg")

    same_samples_img_red = {"LacZ": [],
                            "CG4319": [],
                            "CG12284": []
                            }

    for sample in order:
        for i in order[sample]:
            same_samples_img_red[sample].append([i, red[sample][i]])

    for key in same_samples_img_red:
        if len(same_samples_img_red[key]) == 6:
            get_concat_6(same_samples_img_red[key], key).save(output_path + "merged " + str(key) + "(red).jpg")
        if len(same_samples_img_red[key]) == 3:
            get_concat_3(same_samples_img_red[key], key).save(output_path + "merged " + str(key) + "(red).jpg")
    return True

if __name__ == "__main__":
    images()