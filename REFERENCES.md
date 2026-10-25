# 方法和数据来源

本研究代码实现的是文中明确的改编规则；引用原方法不等于复现其全部研究。

- 生境公式与版本行为：Natural Capital Project. InVEST 3.20.2 Habitat Quality source. https://github.com/natcap/invest/blob/3.20.2/src/natcap/invest/habitat_quality/habitat_quality.py
- 参数语义、威胁权重归一化和质量公式：InVEST Habitat Quality User Guide. https://storage.googleapis.com/releases.naturalcapitalproject.org/invest-userguide/latest/en/habitat_quality.html
- 北京地类生境适宜度、敏感度与最大威胁距离：陈妍、乔飞、江磊，2016，北京大学学报（自然科学版），52(3)，553–562. https://doi.org/10.13209/j.0479-8023.2016.057
- 地类恢复条件系数：Chen et al., 2024, Ecological Indicators 167, 112607. https://doi.org/10.1016/j.ecolind.2024.112607 。本文仅用地类面积比例乘原文系数/100，省略原生态弹性式的多样性乘项，不把系数解释为实际恢复概率。
- 斑块凝聚度：FRAGSTATS C7 Patch Cohesion Index. https://fragstats.org/index.php/fragstats-metrics/patch-based-metrics/aggregation-metrics/c7-patch-cohesion-index 。本文采用/100表达；无目标斑块取0是本研究的明确扩展，缺测不补0。
- CRITIC：Diakoulaki, Mavrotas and Papayannakis, 1995, Computers & Operations Research 22(7), 763–770. https://doi.org/10.1016/0305-0548(94)00059-H 。本文在原有0—1分量上计算标准差与Pearson相关，未做样本极值归一化，仅作权重对照。
- CLCD年度分类：Yang and Huang, Version 1.0.2. https://doi.org/10.5281/zenodo.8176941
- ESA WorldCover 2020 v100：Zanaga et al. https://doi.org/10.5281/zenodo.5571936
- OSM主要道路原始历史提取：Geofabrik China 2020-01-01. https://download.geofabrik.de/asia/china-200101-free.shp.zip 。分布数据遵守其原始署名及ODbL等适用条件。

固定运行环境官方发布页：
https://pypi.org/project/numpy/2.5.4/
https://pypi.org/project/scipy/1.18.1/
https://pypi.org/project/rasterio/1.5.2/
