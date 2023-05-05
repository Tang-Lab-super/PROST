## How to use `PROST`
Before you use `PROST`, you have to make sure that the following two steps are taken place:  
* `PROST_ENV` environment is activated; 
* Add `environment variables` using the following `python code` before using `PROST`：

        import os
        ENVpath = "your path of PROST_ENV"  
        os.environ['R_HOME'] = f'{ENVpath}/lib/R'
        os.environ['R_USER'] = f'{ENVpath}/lib/python3.7/site-packages/rpy2'

**Note: using `conda info -e` can get `your path of PROST_ENV`.**