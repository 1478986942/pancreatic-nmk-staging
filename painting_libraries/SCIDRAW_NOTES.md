# SciDraw SVG 下载 — 用于绘制实验流程图/示意图

## 已下载的素材

| 文件 | 大小 | DOI | 用途建议 |
|------|------|-----|---------|
| `mouse/mouse-profile-f8e60949.svg` | 278 KB | 10.5281/zenodo.3925921 | 小鼠侧面轮廓（大图，适合主图） |
| `mouse/adult-mouse-284dd297.svg` | 44 KB | 10.5281/zenodo.6645614 | 成人小鼠全身（标准比例） |
| `mouse/mouse-head-schema-51ed02be.svg` | 24 KB | 10.5281/zenodo.3925903 | 小鼠头部/脑切片示意 |
| `mouse/mouse-ded16e1d.svg` | 43 KB | 10.5281/zenodo.3925901 | 通用小鼠图标 |
| `mouse/walking-mouse-3c30783c.svg` | 39 KB | 10.5281/zenodo.3925915 | 行走小鼠（适合时间线） |

## 在 R/ggplot2 中使用

```r
library(rsvg)   # 渲染 SVG
library(ggplot2)
library(cowplot)

# 读取 SVG
svg_img <- image_read_svg("painting_libraries/scidraw_assets/mouse/mouse-profile-f8e60949.svg")

# 嵌入 ggplot
p_base <- ggplot(mtcars, aes(wt, mpg)) + geom_point() + theme_classic()
ggdraw(p_base) + draw_image(svg_img, x=0.65, y=0.3, width=0.3, height=0.3)
```

## 在 Python/matplotlib 中使用

```python
from matplotlib.image import imread as mpimg
import matplotlib.pyplot as plt

img = mpimg('painting_libraries/scidraw_assets/mouse/mouse-profile-f8e60949.svg')
fig, ax = plt.subplots(figsize=(8, 5))
ax.imshow(img, extent=[6, 9, 1, 4])   # 放入指定坐标范围
ax.set_xlim(0, 12); ax.set_ylim(0, 8)
plt.savefig('figure_with_mouse.png', dpi=250)
```

## 引用要求

所有 SciDraw 素材使用 CC BY 4.0 协议，必须在论文方法部分引用对应 DOI。
Nature 系列期刊接受 CC BY 4.0 素材，但需在 figure legend 注明：
> "Illustrations in this figure were created using assets from SciDraw (CC BY 4.0)."
