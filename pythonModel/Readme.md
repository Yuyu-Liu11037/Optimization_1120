# Useage

## 1. 安装环境

如果已经单独安装了 CPLEX，请根据官方文档进行配置。以下只提供社区版的安装环境。

创建 Conda 环境，安装依赖：

```
conda install -c ibmdecisionoptimization cplex docplex
```



## 2. 执行带参数命令


无 Callback:
```
python model2.py -n 25 -l n_25_m_0_omega_3.log -o n_25_m_0_omega_3.txt
```

Callback cut:

```
python model2.py -n 25 -l n_25_m_0_omega_3_cut.log -o n_25_m_0_omega_3_cut.txt -c
```


## 3. 其他参数

```
  -c, --cut         To indicate use cut or not
  -n                n
  -m                m
  -e, --omega       omega
  -o, --output      Output file.
  -l, --logout      Output solving log.
  -w, --worker      Number of Threads
  -t, --timelimit   Time Limit in seconds.
```
