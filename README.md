# 汉信码 (Han Xin Code) 纯 Python 编解码器

单文件实现汉信码的**生成**与**解析**，按 **GB/T 21049-2022 / ISO/IEC 20830** 编写，
不依赖任何第三方条码库。

```bash
python hanxin.py selftest                          # 自检
python hanxin.py decode examples/sample-2-generated-text.png
python hanxin.py encode "Hello 汉信码" -o out.png --ec L3
```

## 特点

* **生成**：仅用标准库（自带 PNG 编码器，无需 Pillow）。
* **解析**：需要 `pip install pillow`；自动处理 **90° 整数倍旋转、反色（白码黑底）、
  非整数模块尺寸、缩放与 JPEG 压缩噪声**，并用 RS 纠错恢复受损码字。
* **跨平台**：无 shebang、无 POSIX 专用调用；Windows 控制台（GBK）输出已做 UTF-8 兜底。
  实测 Python 3.11（Windows）/ 3.13（Linux）。
* **完整度**：84 个版本（23×23–189×189）、纠错等级 L1–L4、掩模 0–3；
  解析支持全部数据模式：数字 / 文本 / 常用汉字一区二区 / GB 18030 双·四字节 / 字节 /
  ECI / Unicode（差分 UTF-8）/ GS1 / URI。

## 安装

```bash
git clone https://github.com/CPPU-ISA/HanXinCodec.git
cd HanXinCodec
pip install pillow          # 仅"解析图片"需要；只生成可不装
```

## 命令行

```bash
python hanxin.py selftest                                  # 自检：矩阵往返 + 图片往返
python hanxin.py decode code.png                           # 解析（自动判版本/旋转/反色）
python hanxin.py decode code.png -j                        # 输出 JSON（分段、RS 块、符号标识符）
python hanxin.py decode code.png --ascii                   # 同时打印模块矩阵
python hanxin.py decode matrix.txt                         # 也可解析 0/1 或 '#'/'.' 文本矩阵

python hanxin.py encode "文本" -o out.png                  # 自动选最小版本 + 最优掩模
python hanxin.py encode "文本" -o out.png --ec L4 --version 6 --mask 2 --scale 12
```

## 作为库调用

```python
import hanxin

# 生成
result = hanxin.encode("汉信码测试", error_correction="L2")
print(result["version"], result["dimension"], result["mask"])       # 1 23 1
png = hanxin.render_png(result["modules"], result["dimension"], scale=10)
open("out.png", "wb").write(png)

# 解析
info = hanxin.decode_image("out.png")
print(info["text"], info["error_level"], info["corrected_errors"])

# 直接对矩阵解码（list[list[int]] 或扁平 list 均可）
info = hanxin.decode_matrix(result["modules"], result["dimension"])
```

主要 API：`encode()`、`decode_matrix()`、`decode_image()`、`image_to_matrix()`、
`render_png()`、`rs_encode()`、`rs_decode()`、`rs_structure(version, level)`。

## 示例图

| 文件 | 解析结果 | 说明 |
| --- | --- | --- |
| `examples/sample-1-url.png` | `https://mp.weixin.qq.com/s/bLn1uNnn_0e3mI1BpMMgww` | 版本 7（35×35），模块 13.2px **非整数**，实测 0 纠错 |
| `examples/sample-2-generated-text.png` | `Hello Han Xin Code 2026` | 本脚本生成，文本 + 数字分段 |
| `examples/sample-3-generated-chinese.png` | `汉信码 Han Xin Code 2026` | 本脚本生成，常用汉字 + 文本 + 数字混合分段 |

## 交叉验证

| 验证项 | 结果 |
| --- | --- |
| 码块结构 / 对齐参数 vs Zint（BSD-3-Clause） | 84 版本 × 4 等级 = 336 组全部一致 |
| Zint 生成 → 本实现解析（版本 1–84 全扫） | 84/84 文本、版本、等级、掩模全对 |
| 本实现生成 → 独立 JS 参考实现解析 | 全部正确，0 纠错 |
| 旋转 / 反色 / 非整数缩放 / 150px 小图 / JPEG q30 | 全部正确，0 纠错 |

## 限制

* 图片前端假设符号**轴对齐、静区干净**；**不支持任意角度倾斜与透视畸变**的照片
  （这类场景需要角点检测 + 单应性变换）。
* 生成端未实现 URI / GS1 / ECI / Unicode 的**编码**（解析端全部支持），
  此类内容会退化为字节模式，体积略大但解码结果一致。

## 许可与来源

**MIT**（见 `LICENSE`）。实现为独立编写的代码，未包含任何第三方源码；
文件内的参数表是标准的**规范性数据**（互操作所必需），已与 Zint 逐项核对一致；
仓库不包含标准正文（GB/T 与 ISO 标准文本受版权保护，请勿随代码分发）。
