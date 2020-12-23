import sys
import datetime
class Logger():
    def __init__(self,filename,stream=sys.stdout):
        # filename指明是train，eval还是benchmark
        self._filename='./logs/'+datetime.datetime.now().strftime('%Y_%m_%d_%H_%M_%S')+"_"+filename+"_log.txt"
        self._terminal=stream
        self._log=open(self._filename,"a")
    def write(self,message):
        # 实现终端输出和log文件保存
        self._terminal.write(message)
        self._log.write(message)
    def flush(self):
        pass