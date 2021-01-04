import numpy as np
from scipy.special import lambertw
import matplotlib.pyplot as plt
import collections
import json
import os
class Cnoise():
    def __init__(self,hnum):
        self._samples=hnum
    
    def noisepower(self,x,snr):
        # sigpower=np.sum(np.abs(x)**2,axis=(0,1))/(x.shape[0]*x.shape[1])
        
        sigpower=np.sum(x**2)/(x.shape[0]*x.shape[1])
        snr=10**(snr/10)
        npower=sigpower/snr
        return npower
    def gaussion(self,x,npower):
        row,col=x.shape
        # snr = 10 ** (snr / 10.0)
        # xpower = np.sum(x ** 2) / len(x)
        # npower = xpower / snr
        noise_real = np.random.randn(row,col) * np.sqrt(npower)
        noise_imag = np.random.randn(row,col) * np.sqrt(npower)
        noise=noise_real+1j*noise_imag
        return noise
class gaussionnoise(Cnoise):
    def __init__(self,hnum):
        super().__init__(hnum)
        print("gaussion noise generate")
    def getawgn(self,x,snr):
        npower=self.noisepower(x,snr)
        noise=self.gaussion(x,npower)
        return noise
class itsnoise(Cnoise):
    def __init__(self,config,hnum):
        super().__init__(hnum)
        print("its noise generate")
        self._daqizaoshengfile0="data/noise/daqizaosheng0.json"
        self._daqizaoshengfile1="data/noise/daqizaosheng1.json"
        self._count=0
        self._Nr=config["Nr"]
        self._timeseq=np.arange(0,hnum*config["sample_interval"],config["sample_interval"]).reshape((-1,1))
        self._noiseconfig={
            "carrierfreq":13.86e6,
            "lpbw":config["bandwidth"],
            "theta_zg":2.4,
            "gamma_zg":0.15,
            "Ni":40,
            "theta_mg":1.2,
            "gamma_mg":1e-8,
            "Nj":200,
            "mg_max":2.0e-5,
            "mg_section":(450e-6,550e-6),
            "mg_interval":4e-6,
            "theta_dz":5,
            "gamma_dz":1,
            "z1_tf":57.43,
            "z2_tf":32.23,
            "z3_tf":12.68,
            "z1_aj":18.62,
            "z2_aj":16.62,
            "z3_aj":1.49,
            "sample_interval":config["sample_interval"],
        }
        print("channel noise config",self._noiseconfig)
        self._allsnrdz=[]
        # if not os.path.exists(self._daqizaoshengfile0):
        #     self.daqizaoshengsave(self._daqizaoshengfile0)
        # if not os.path.exists(self._daqizaoshengfile1):
        #     self.daqizaoshengsave(self._daqizaoshengfile1)
        # with open(self._daqizaoshengfile0,"r",encoding='utf=8') as f:
        #     self._allsnrdz.append(json.load(f))
        # with open(self._daqizaoshengfile1,"r",encoding='utf=8') as f:
        #     self._allsnrdz.append(json.load(f))
    def hall1(self,li,theta,gamma):
        res=gamma*np.sqrt(np.power(1-li,1/(1-theta))-1)
        return res
    def hall2(self,li,theta,gamma,mgmax):
        res=gamma*np.sqrt(np.power((li*(np.power(mgmax**2/gamma**2+1,(1-theta)/2)-1)+1),2/(1-theta))-1)
        return res
    def hall3(self,li,z1,z2,z3):
        res=-1*z1/(z2*z3)-1/z3*np.log(1-li)+1/z2*lambertw(z1/z3*np.exp((z1+z2*np.log(1-li))/z3))
        return res
    def zhaidaiganrao(self):
        # gamma_zg=np.sqrt(0.5*(-npower+np.sqrt(npower**2+4*(self._noiseconfig["theta_zg"]-1)/(self._noiseconfig["theta_zg"]+1))))
        # print("gamma_zg=",gamma_zg,"npower=",npower)
        gamma_zg=self._noiseconfig["gamma_zg"]
        zg=np.zeros((self._samples,1))+1j*np.zeros((self._samples,1))
       
        for index in range(self._samples):
            for i in range(0,self._noiseconfig["Ni"]):
                distribution=np.random.uniform(0,1)
                zg_maganitude=self.hall1(distribution,self._noiseconfig["theta_zg"],gamma_zg)
                omega=np.random.uniform(-1*self._noiseconfig["lpbw"],self._noiseconfig["lpbw"])
                fai=np.random.uniform(0,np.pi*2)
                zg[index,0]+=zg_maganitude*np.exp(-1j*(omega*self._timeseq[index]+fai))
            if np.abs(np.real(zg[index,0]))>50 or np.abs(np.imag(zg[index,0]))>50:
                zg[index,0]=0
        return zg
    
    def maichongganrao(self):
        mg=np.zeros((self._samples,1))+1j*np.zeros((self._samples,1))
        # gamma_mg=np.sqrt(0.5*(-npower+np.sqrt(npower**2+4*(self._noiseconfig["theta_mg"]-1)/(self._noiseconfig["theta_mg"]+1))))
        # print("gamma_mg=",gamma_mg,"npower=",npower)
        gamma_mg=self._noiseconfig["gamma_mg"]
        mg_tj=[0]
        while len(mg_tj)<=self._noiseconfig["Nj"]:
            t=np.random.uniform(self._noiseconfig["mg_section"][0],self._noiseconfig["mg_section"][1])
            t+=mg_tj[-1]
            mg_tj.append(np.random.uniform(t,t+self._noiseconfig["mg_interval"]))
        mg_tj.pop(0)
        for index in range(self._samples):
            for i in range(self._noiseconfig["Nj"]):
                distribution=np.random.uniform(0,1)
                mg_maganitude=self.hall2(distribution,self._noiseconfig["theta_mg"],gamma_mg,self._noiseconfig["mg_max"])
                mg[index,0]+=mg_maganitude*np.sin(2*np.pi*self._noiseconfig["lpbw"]*(self._timeseq[index]-mg_tj[i]))/(self._timeseq[index]-mg_tj[i])*np.exp(1j*self._noiseconfig["carrierfreq"]*mg_tj[i])
        return mg
    def daqizaosheng(self,snr,index):
        dz_maganitude=self._allsnrdz[index][str(snr)]
        dz=np.zeros((self._samples,1))+1j*np.zeros((self._samples,1))
        for index in range(self._samples):
            fai=np.random.uniform(0,1)
            dz[index,0]=dz_maganitude[index]*np.exp(1j*fai)
        return dz
    def daqizaoshengsave(self,filename):
        sigpower=1
        snrlist=[i for i in range(-10,31,2)]
        dzm_dic={}
        # 所有snr具有相同的突发期和安静期
        distribution=np.random.uniform(0,1)
        T_tf=self.hall3(distribution,self._noiseconfig["z1_tf"],self._noiseconfig["z2_tf"],self._noiseconfig["z3_tf"])
        distribution=np.random.uniform(0,1)
        T_aj=self.hall3(distribution,self._noiseconfig["z1_aj"],self._noiseconfig["z2_aj"],self._noiseconfig["z3_aj"])
        num_tf=int(T_tf/self._noiseconfig["sample_interval"])
        num_aj=int(T_aj/self._noiseconfig["sample_interval"])

        for snr in snrlist:
            snrlinear=10**(snr/10)
            npower=sigpower/snrlinear
            gamma_dz=np.sqrt(0.5*(-npower+np.sqrt(npower**2+4*(self._noiseconfig["theta_dz"]-1)/(self._noiseconfig["theta_dz"]+1))))
            pshreshold=float(gamma_dz**2*(np.power(1+T_aj/T_tf,2/(self._noiseconfig["theta_dz"]-1))-1))
            print("snr=",snr,"gamma_dz=",gamma_dz,"shreshold=",pshreshold)
            dz_maganitude=[]
            def genedz(state):
                distribution=np.random.uniform(0,1)
                res=self.hall1(distribution,self._noiseconfig["theta_dz"],gamma_dz)
                while (state=="tf" and res<=pshreshold) or (state=="aj" and res>=pshreshold):
                    distribution=np.random.uniform(0,1)
                    res=self.hall1(distribution,self._noiseconfig["theta_dz"],gamma_dz)
                return res

            while len(dz_maganitude)<self._samples:
                count=0
                while len(dz_maganitude)<self._samples and count<num_tf:
                    count+=1
                    dz_maganitude.append(genedz("tf"))
                count=0
                while len(dz_maganitude)<self._samples and count<num_aj:
                    count+=1
                    dz_maganitude.append(genedz("aj"))
            dzm_dic[str(snr)]=dz_maganitude
        
        with open(filename,"w",encoding='utf=8') as f:
            json.dump(dzm_dic,f,ensure_ascii=False,indent=2)
        

    def getits(self):
        # npower=self.noisepower(x,snr)
        # print("***get its noise***")
        zg=self.zhaidaiganrao()
        mg=self.maichongganrao()
        # dz=self.daqizaosheng(snr,self._count)
        self._count=(self._count+1)%self._Nr #用来区分是第一个接收机
        # plt.figure(1)
        # plt.plot(self._timeseq,np.imag(zg))
        # plt.figure(2)
        # plt.plot(self._timeseq,np.imag(mg))
        # # plt.figure(3)
        # # plt.plot(self._timeseq,np.imag(dz))
        # plt.show()
        return zg+mg
    def getawgn(self,x,snr):
        npower=self.noisepower(x,snr)
        noise=self.gaussion(x,npower)
        return noise
    
if __name__ == "__main__":
    sps=1
    bandwidth=250e3
    roll_off=0.25
    baud_rate=bandwidth/(1+roll_off)
    sample_rate=baud_rate/sps
    sample_interval=1/sample_rate
    params={
        "sps":sps,
        "bandwidth":bandwidth,
        "roll_off":roll_off,
        "baud_rate":baud_rate,
        "sample_rate":sample_rate,
        "sample_interval":sample_interval,
        "Nr":2
    }
    hnum=10000
    noise=itsnoise(params,hnum)
    x=np.random.randint(0,2,(hnum,1))
    xpower = np.sum(x ** 2) / len(x)
    for i in range(hnum):
        if x[i,0]==0:
            x[i,0]=-1
    snrlist=[i for i in range(-10,31,2)]
    for snr in snrlist:
        a=noise.getits()
        x=noise.getawgn(x,snr)


            




    
