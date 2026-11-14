# 本论文 GitHub 复现材料

本仓库提供巴西与墨西哥短期高温及结石相关住院记录研究的复现材料。仓库：https://github.com/FENG1567/latam-heat-stone-admissions；版本：v1.0.0。代码可用性声明见CODE_AVAILABILITY.md。

## 可以直接完成的事

1. 读取校验通过的最终冻结结果，重新生成4幅正文图、补图S1、4张主表及S1–S16补充表。
2. 获取论文87个源CSV的科学名称映射、字段表、冻结结果与最终稿的版本关系。
3. 无损还原随附的地区/州日与临床分层汇总输入，供统计重拟合使用。大文件按字节分片；还原后SHA256与校正后的服务器原文件完全一致，不是重新聚合或重算。
4. 从保留的分析源代码创建独立运行目录，检查依赖和所需输入，再按模块执行。

## 不应误解的地方

绘图复现不是统计重拟合；统计重拟合也不是从所有原始行政文件和气象格点重新提取。后者还涉及第三方原始文件、历史检索和外部地理文件，不能把这些未重新验证的环节说成已一键跑通。本次没有重跑全部全国模型，也没有把数十GB历史拟合缓存上传。具体已验证范围、手动前置项见docs/REPRODUCIBILITY.md和模块清单。

## 最简命令

```sh
python -m pip install -r environment/requirements-publication.txt
python scripts/verify_repo.py
python code/publication/run_publication.py --frozen data/frozen/frozen_numbers.md --outdir build/publication
python scripts/analysis.py prepare
python scripts/analysis.py list
python scripts/analysis.py preflight
```

完整操作见README.md。首次打开请先看docs/UPLOAD_GITHUB_中文.md。自编代码MIT、自有汇总结果CC BY 4.0已获作者确认；第三方原始和衍生输入保留来源条款。

墨西哥最终队列runner在models和smoke阶段也会先从所附来源汇总重建日历，并可能改写工作副本的压缩字节；仓库中的原始分片不受影响。模型运行后如需再次按原SHA还原，请指定新的独立workspace，不要覆盖已改变的工作目录。
