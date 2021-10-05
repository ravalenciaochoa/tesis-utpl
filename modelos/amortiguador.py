import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt
from openseespy.postprocessing.Get_Rendering import plot_model


class AmortiguadorViscoso():

    
    #Variables globales : CONSTANTES

    PERIODO=0.7
    KD=25
    CD=20.7452
    AD=0.35
    H=3
    PESO=10000

    
    #RECODERS
    FILE_DESPLAZAMIENTO_N1 = "Despla_nodo_1.out"    # DESPLAZAMIENTO_N1
    FILE_DESPLAZAMIENTO_N2 = "Despla_nodo_2.out"    # DESPLAZAMIENTO_N2
    FILE_REAC_BASE = "Reacc_Base.out"               # REACCION_BASE      
    FILE_FUERZAS_COL = "Fuerzas_Col.out"            # FUERZAS_COLUMNA
    FILE_FUERZAS_COL_SEC = "Fuerzas_Col_sec.out"    # FUERZAS_COLUMNA SEC
    FILE_DEFORMACION_COL_SEC = "Defor_col_sec.out"  # DEFORMACION



    def __init__(self, nombre='BasicBuilder', numero_dimensiones=2, numero_dots=3,altura_columna=0.6, base_columna=0.6, numero_modos=2,data_file="amortiguador.txt",calculo='*'): 
        
        self.numero_dimensiones = numero_dimensiones
        self.numero_dots = numero_dots
        self.altura_columna = altura_columna
        self.base_columna = base_columna
        self.nombre = nombre
        self.numero_modos = numero_modos
        self.data_file = data_file
        self.calculo = calculo

        #Creación del modelo 
        ops.wipe()
        ops.model(self.nombre, "-ndm", self.numero_dimensiones, "-ndf", self.numero_dots)

    def construir_modelo(self):

        #Definir los nodos de la estructura
        ops.node(1, 0, 0)
        ops.node(2, 0, self.H)

        #Restricciones    
        ops.fix(1, 1, 1, 1)
        ops.fix(2, 0, 1, 0)

        #Definir la masa en la cabeza de la columna
        ops.mass(2, self.PESO, 1e-8, 1e-8)

        #Definir las propiedades de la columna
        #K=(np.pi/self.PERIODO)**2
        E=20000000
        Ic=(self.H**3/24*E)

        area = self.altura_columna*self.base_columna

        #Definir el material viscoso
        ops.uniaxialMaterial("ViscousDamper", 1, self.KD, self.CD, self.AD)

        #Definir la transformacion de coordenadas
        transformacion_coordenadas=1
        ops.geomTransf("Linear", transformacion_coordenadas)

        #Se definen los elementos.
        ops.element("elasticBeamColumn", 1, 1, 2, area, E, Ic, transformacion_coordenadas)

        # Amortiguador
        #element twoNodeLink $eleTag $iNode $jNode -mat $matTags -dir $dirs
        ops.element("twoNodeLink", 2, 1, 2, "-mat", 1, "-dir", 1)

        #Definir los recorders
        ops.recorder("Node", "-file", self.FILE_DESPLAZAMIENTO_N2, "-time", "-node", 2, "-dof", 1, 2, 3, "disp")
        ops.recorder("Node", "-file", self.FILE_DESPLAZAMIENTO_N1, "-time", "-node", 1, "-dof", 1, 2, 3, "disp")
        ops.recorder("Node", "-file", self.FILE_REAC_BASE, "-time", "-node", 1, "-dof", 1, 2, 3, "reaction")
        ops.recorder("Element", "-file", self.FILE_FUERZAS_COL, "-time", "-ele", 1, "localForce")
        ops.recorder("Element", "-file", self.FILE_FUERZAS_COL_SEC, "-time", "-ele", 1, "deformations")
        ops.recorder("Element", "-file", self.FILE_DEFORMACION_COL_SEC, "-time", "-ele", 1, "-dof", 1, "force")

        #Analisis modal
        #numero_modos=20 #int(input("Ingrese el numero de modos a calcular: "))
        
        omega_sq = ops.eigen("-fullGenLapack", self.numero_modos)
        omega_sq = np.array(omega_sq)

        omega = omega_sq**0.5
        periodo =(2*np.pi)/omega

        for i in range(self.numero_modos):
            print("Modo", i+1, "T= ", periodo[i], "s")

        #Definir análisis dinamico
        factor=0.001
        dt=0.001
        n_puntos=15000
        
        var_dir=1
        ops.timeSeries("Path", 2, "-dt", dt, "-filePath", self.data_file, "-factor", factor)
        ops.pattern("UniformExcitation", 2, var_dir, "-accel", 2)

        #Parametros de analisis dinamico
        ops.constraints("Transformation")
        ops.numberer("RCM")
        ops.system("UmfPack")
        ops.test("EnergyIncr", 1.0e-8, 10)
        ops.algorithm("KrylovNewton")
        ops.integrator("Newmark", 0.5, 0.25)
        ops.analysis("Transient")
        ops.analyze(n_puntos, dt)
        ops.remove("recorders")


        if self.calculo == 'DESPLAZAMIENTO_N2' :    
            array_mek=np.loadtxt(self.FILE_DESPLAZAMIENTO_N2)            
            self.construye_figura(array_mek,serie=[1,2,3], recoder=self.calculo)                            
        elif self.calculo == 'DESPLAZAMIENTO_N1':
            array_mek=np.loadtxt(self.FILE_DESPLAZAMIENTO_N1)                
            self.construye_figura(array_mek,serie=[1,2,3], recoder=self.calculo)                            
        elif self.calculo == 'REACCION_BASE':
            array_mek=np.loadtxt(self.FILE_REAC_BASE)
            self.construye_figura(array_mek,serie=[1,2,3], recoder=self.calculo)                            
        elif self.calculo == 'FUERZAS_COLUMNA':
            array_mek=np.loadtxt(self.FILE_FUERZAS_COL)
            self.construye_figura(array_mek,serie=[1], recoder=self.calculo)                            
        elif self.calculo == 'FUERZAS_COLUMNA_SEC':
            array_mek=np.loadtxt(self.FILE_FUERZAS_COL_SEC)
            self.construye_figura(array_mek,serie=[1], recoder=self.calculo)                            
        elif self.calculo == 'DEFORMACION_COL_SEC':
            array_mek=np.loadtxt(self.FILE_DEFORMACION_COL_SEC)            
            self.construye_figura(array_mek,serie=[1], recoder=self.calculo)                            
        else:
            datos_1=np.loadtxt(self.FILE_DESPLAZAMIENTO_N2)                                    
            datos_2=np.loadtxt(self.FILE_DESPLAZAMIENTO_N1)                            
            datos_3=np.loadtxt(self.FILE_REAC_BASE)            
            datos_4=np.loadtxt(self.FILE_FUERZAS_COL)            
            serie_total = [[ datos_1, [1,2,3],'DESPLAZAMIENTO_N2' ], [ datos_2, [1,2,3],'DESPLAZAMIENTO_N1' ], [ datos_3, [1,2,3],'REACCION_BASE' ], [ datos_4, [1],'FUERZAS_COLUMNA' ] ]     

            if serie_total:
                self.construye_figura(serie_total=serie_total, recoder=self.calculo)
                        

    def construye_figura(self, array_mek=None, serie=[], recoder=None, serie_total=[], resultados=[]):
            
        """Construye una figura basado en un array de datos, una serie de gráficos        
        """

        if recoder != '*':

            nombre_figura = "Gráfico # "
            titulo = None

            for i in serie:

                if recoder != None:
                    titulo = recoder+' : '+nombre_figura+str(i)
                    plt.figure(titulo)
                else:
                    titulo = nombre_figura + str(i)
                    plt.figure(titulo)
                
                plt.plot(array_mek[:,0], array_mek[:,i])
                
                plt.ylabel('Desplazamiento Horizontal (m)')
                plt.xlabel('Tiempo (s)')                       
                plt.show()

        else:
            from PIL import Image
            import matplotlib.image as mpimg
            import os
            from pathlib import Path
            
            BASE_DIR = Path(__file__).resolve().parent.parent
            IMAGEN_COLUMNA = os.path.join(os.path.join(BASE_DIR,'modelos'), 'columna_con_carga_ciclica.jpeg')
            pil_img = Image.open(IMAGEN_COLUMNA)
                                                
            fig = plt.figure()


            FILAS = 3
            COLUMNAS = 4

            ax1 = plt.subplot2grid((FILAS, COLUMNAS), (0, 0), rowspan=2)
            ax1.imshow(pil_img) 
            ax1.axes.xaxis.set_visible(False)
            ax1.axes.yaxis.set_visible(False)    

            ax2 = plt.subplot2grid((FILAS, COLUMNAS), (2, 0))

            ax3 = plt.subplot2grid((FILAS, COLUMNAS), (0, 1))
            ax4 = plt.subplot2grid((FILAS, COLUMNAS), (1, 1))
            ax5 = plt.subplot2grid((FILAS, COLUMNAS), (2, 1))

            ax6 = plt.subplot2grid((FILAS, COLUMNAS), (0, 2))
            ax7 = plt.subplot2grid((FILAS, COLUMNAS), (1, 2))
            ax8 = plt.subplot2grid((FILAS, COLUMNAS), (2, 2))

            ax9 = plt.subplot2grid((FILAS, COLUMNAS), (0, 3))
            ax10 = plt.subplot2grid((FILAS, COLUMNAS), (1, 3))
            ax11 = plt.subplot2grid((FILAS, COLUMNAS), (2, 3))
            
            
            for columna, elemento in enumerate(serie_total):                
                
                array_mek = elemento[0]
                serie = elemento[1]
                calculo_var = elemento[2]
            
                if calculo_var == 'DESPLAZAMIENTO_N2':  
                    
                    ax3.plot(array_mek[:,0], array_mek[:,serie[0]],'tab:red')
                    ax3.set_title(f" {calculo_var}\nGRAFICO : # {serie[0]}", fontsize=8)
                    ax3.set_xlabel('Tiempo (seg)')
                    ax3.set_ylabel('Desplazamiento H. (m)')
                                        
                    ax4.plot(array_mek[:,0], array_mek[:,serie[1]],'tab:red')
                    ax4.set_title(f" {calculo_var}\nGRAFICO : # {serie[1]}", fontsize=8)
                    ax4.set_xlabel('Tiempo (seg)')
                    ax4.set_ylabel('Desplazamiento H. (m)')
                    
                    ax5.plot(array_mek[:,0], array_mek[:,serie[2]],'tab:red')
                    ax5.set_title(f" {calculo_var}\nGRAFICO : # {serie[2]}", fontsize=8)
                    ax5.set_xlabel('Tiempo (seg)')
                    ax5.set_ylabel('Desplazamiento H. (m)')

                elif calculo_var == 'DESPLAZAMIENTO_N1':
                                        
                    ax6.plot(array_mek[:,0], array_mek[:,serie[0]],'tab:orange')
                    ax6.set_title(f" {calculo_var}\nGRAFICO : # {serie[0]}", fontsize=8)
                    ax6.set_xlabel('Tiempo (seg)')
                    ax6.set_ylabel('Desplazamiento H. (m)')
                                        
                    ax7.plot(array_mek[:,0], array_mek[:,serie[1]],'tab:orange')
                    ax7.set_title(f" {calculo_var}\nGRAFICO : # {serie[1]}", fontsize=8)
                    ax7.set_xlabel('Tiempo (seg)')
                    ax7.set_ylabel('Desplazamiento H. (m)')
                    
                    ax8.plot(array_mek[:,0], array_mek[:,serie[2]],'tab:orange')
                    ax8.set_title(f" {calculo_var}\nGRAFICO : # {serie[2]}", fontsize=8)
                    ax8.set_xlabel('Tiempo (seg)')
                    ax8.set_ylabel('Desplazamiento H. (m)')

                elif calculo_var == 'REACCION_BASE':
                                        
                    ax9.plot(array_mek[:,0], array_mek[:,serie[0]],'tab:blue')
                    ax9.set_title(f" {calculo_var}\nGRAFICO : # {serie[0]}", fontsize=8)
                    ax9.set_xlabel('Tiempo (seg)')
                    ax9.set_ylabel('Desplazamiento H. (m)')
                                        
                    ax10.plot(array_mek[:,0], array_mek[:,serie[1]],'tab:blue')
                    ax10.set_title(f" {calculo_var}\nGRAFICO : # {serie[1]}", fontsize=8)
                    ax10.set_xlabel('Tiempo (seg)')
                    ax10.set_ylabel('Desplazamiento H. (m)')
                    
                    ax11.plot(array_mek[:,0], array_mek[:,serie[2]],'tab:blue')
                    ax11.set_title(f" {calculo_var}\nGRAFICO : # {serie[2]}", fontsize=8)
                    ax11.set_xlabel('Tiempo (seg)')
                    ax11.set_ylabel('Desplazamiento H. (m)')

                else :

                    ax2.plot(array_mek[:,0], array_mek[:,serie[0]],'tab:green')  
                    ax2.set_title(f" {calculo_var}\nGRAFICO : # {serie[0]}", fontsize=8)
                    ax2.set_xlabel('Tiempo (seg)')
                    ax2.set_ylabel('Desplazamiento H. (m)')
                
            fig.suptitle("MODELOS MATEMÁTICOS: AMORTIGUADOR VISCOSO ")                             
            fig.subplots_adjust(left=0.067, bottom=0.107, right=0.971, top=0.88, wspace=0.562, hspace=1)

            #Para fijar un titulo a la ventana principal
            mng = plt.get_current_fig_manager()
            mng.set_window_title("MODELOS MATEMÁTICOS UTPL")
            mng.resize(*mng.window.maxsize())                            
            plt.show() 