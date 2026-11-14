# GitHub仓库与版本

公开仓库：https://github.com/FENG1567/latam-heat-stone-admissions

版本：v1.0.0。自编代码使用MIT；自有汇总结果和图件使用CC BY 4.0；第三方输入保留提供方原有条款。

获取版本：

```sh
git clone --branch v1.0.0 https://github.com/FENG1567/latam-heat-stone-admissions.git
cd latam-heat-stone-admissions
python -m pip install -r environment/requirements-publication.txt
python scripts/verify_repo.py
python code/publication/run_publication.py --frozen data/frozen/frozen_numbers.md --outdir build/publication
```

统计重拟合使用analysis入口和模块清单，具体依赖与外部输入见README.md和docs/REPRODUCIBILITY.md。图表复现选择冻结估计，不拟合模型。源码、汇总输入和源CSV均保留目录结构。

大输入使用二进制分片，运行prepare_inputs.py还原；请勿手工改写分片或冻结表。保留.gitattributes以避免跨平台换行转换破坏校验值。build、凭据、患者级原始记录和拟合缓存不属于发布内容。
