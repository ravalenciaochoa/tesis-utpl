import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt
from openseespy.postprocessing.Get_Rendering import plot_model

class ColumnaAislador():

    FILE_DESPLAZAMIENTO_NODO2 = "DESPLAZAMIENTO_NODO2.out"          # DESPLAZAMIENTO N2
    FILE_DESPLAZAMIENTO_NODO1 = "DESPLAZAMIENTO_NODO1.out"          # DESPLAZAMIENTO N1
    FILE_REAC_BASE = "REAC_BASE.out"                                # REACCION_BASE      
    FILE_FUERZAS_COL = "FUERZAS_COL.out"                            # FUERZAS COL
    FILE_FUERZAS_COL_N1 = "FUERZAS_COL_N1.out"                      # FUERZAS COL_N1
    FILE_DEFORMACION_COL = "DEFORMACION_COL.out"                    # DEFORMACION_COL
    FILE_FORCE_COLUMNA = "Data/ForceColSecnumIntgrPts.out"          # FORCE_COLUMNA      
    FILE_DEFORMACION_COLUMNA = "Data/DefoColSec$numIntgrPts.out"    # DEFORMACION_COLUMNA


    def __init__(self,nombre='BasicBuilder',numero_dimensiones=2, numero_dots=3, altura_columna=2, fuerza_aplicada=15000 , resistencia_fc=21000 ,dim_travs_colum_altura=0.5, dim_trav_colum_base=0.5, calculo='*', data_file='datos.txt'):
        
        self.numero_dimensiones = numero_dimensiones
        self.numero_dots = numero_dots
        self.nombre = nombre
        self.dim_travs_colum_altura = dim_travs_colum_altura
        self.dim_trav_colum_base = dim_trav_colum_base
        self.calculo = calculo
        self.altura_columna = altura_columna
        self.fuerza_aplicada = fuerza_aplicada
        self.resistencia_fc = resistencia_fc
        self.data_file = data_file

        ops.wipe()
        #Creación del modelo 
        ops.model(self.nombre, "-ndm", self.numero_dimensiones, "-ndf", self.numero_dots)

    
    def construir_modelo(self):

        #Datos de la estructura                                        
        pes_esp_hor=23.56
        g=9.81
        E=20000000
        col_sec_tag=1
        
        peso_ele=pes_esp_hor*self.dim_travs_colum_altura*self.dim_trav_colum_base*self.altura_columna
        peso_estru=peso_ele+self.fuerza_aplicada
        
        masa=peso_estru/g
        masa_1=masa-323.896
        area=self.dim_travs_colum_altura*self.dim_trav_colum_base
        I=(1/12)*self.dim_trav_colum_base*self.dim_travs_colum_altura**3
    

        #Definir los nodos de la estructura
        ops.node(1, 0, 0)
        ops.node(2, 0, self.altura_columna)

        #Definir los nodos del aislador
        ops.node(3, 0, 0)

        #Definir las restricciones
        ops.fix(1, 0, 1, 1) # restricción para el pie de columna
        ops.fix(3, 1, 1, 1) # restricción para el aislador

        # Definir las secciones elasticas de la columna
        ops.section("Elastic", col_sec_tag, E, area, I)

        # Definir las secciones elasticas del aislador
        ops.uniaxialMaterial("Elastic", 4, 20, 0.1)
        ops.uniaxialMaterial("Elastic", 5, 21000000)

        #Definir la transformacion de coordenadas
        ops.geomTransf("Linear", 1, 0, 0, 1)

        # Definir elementos de columna
        ops.element("elasticBeamColumn", 6, 1, 2, area, E, I, 1)

        # Definir elementos del aislador
        ops.element("zeroLength", 7, 3, 1, "-mat", 4, 5, "-dir", 1, 2)

        #Definir la masa de la estructura
        ops.mass(2, masa_1, 0, 0)

        # Definir los recorders

        ops.recorder("Node", "-file", self.FILE_DESPLAZAMIENTO_NODO2, "-time", "-node", 2, "-dof", 1,2,3, "disp")
        ops.recorder("Node", "-file", self.FILE_DESPLAZAMIENTO_NODO1, "-time", "-node", 1, "-dof",1,2,3 , "disp")
        ops.recorder("Node", "-file", self.FILE_REAC_BASE, "-time", "-node", 1, "-dof", 1, 2, 3, "reaction")
        ops.recorder("Element", "-file", self.FILE_FUERZAS_COL, "-time", "-ele", 6, "localForce")

        #Definir las actuantes sobre la columna
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.eleLoad("-ele", 6, "-type", "-beamUniform", 0, 0, -peso_ele)
        ops.load(2, 0, masa, 0)

        #Parametros de analisis estatico
        ops.constraints("Plain")
        ops.numberer("RCM")
        ops.system("BandGeneral")
        ops.test("NormDispIncr", 1.0e-8, 6)
        ops.algorithm("Newton")

        numero_pasos=10
        ops.integrator("LoadControl", 1/numero_pasos)
        ops.analysis("Static")
        ok=ops.analyze(numero_pasos)

        if ok==0:
            print("CARGA AXIAL APLICADA")
        else:
            print("ERROR AL APLICAR LA CARGA AXIAL")
            
        ops.loadConst("-time", 0.0)

        #Definir el amortiguamiento
        num_modos=1
        omega_sq=ops.eigen("-fullGenLapack", num_modos)
        omega_sq=np.array(omega_sq)
        omega=omega_sq**0.5

        periodo=(2*np.pi)/omega
        for i in range(num_modos):
            print("Modo", i+1, "T= ", periodo[i], "s")
            
        #Definir análisis dinamico
        factor=0.01
        dt=0.01

        var_dir=1
        ops.timeSeries("Path", 2, "-dt", dt, "-filePath", self.data_file, "-factor", factor)
        ops.pattern("UniformExcitation", 2, var_dir, "-accel", 2)

        #Parametros de analisis dinamico
        ops.constraints("Transformation")
        ops.numberer("RCM")
        ops.system("UmfPack")
        ops.test("EnergyIncr", 1.0e-6, 6, 0)
        ops.algorithm("ModifiedNewton")
        ops.integrator("Newmark", 0.5, 0.25)
        ops.analysis("Transient")
        numero_puntos=50000
        ops.analyze(numero_puntos, dt)
        ops.remove("recorders")



        if self.calculo == 'REACCION_BASE' :                
            array_mek=np.loadtxt(self.FILE_REAC_BASE)
            self.construye_figura(array_mek,serie=[1,2,3], recoder=self.calculo)                            

        elif self.calculo == 'DESPLAZAMIENTO_NODO2':            
            array_mek=np.loadtxt(self.FILE_DESPLAZAMIENTO_NODO2)
            self.construye_figura(array_mek,serie=[1,2,3], recoder=self.calculo)

        elif self.calculo == 'DESPLAZAMIENTO_NODO1':            
            array_mek=np.loadtxt(self.FILE_DESPLAZAMIENTO_NODO1)
            self.construye_figura(array_mek,serie=[1,2,3], recoder=self.calculo)
        
        else:            
            datos_1=np.loadtxt(self.FILE_REAC_BASE)
            #self.construye_figura(array_mek,serie=[1,2,3], recoder='REACCION_BASE')

            datos_2=np.loadtxt(self.FILE_DESPLAZAMIENTO_NODO2)
            #self.construye_figura(array_mek,serie=[1,2,3], recoder='DESPLAZAMIENTO_NODO2')

            datos_3=np.loadtxt(self.FILE_DESPLAZAMIENTO_NODO1)
            #self.construye_figura(array_mek,serie=[1,2,3], recoder='DESPLAZAMIENTO_NODO1')

            serie_total = [ [ datos_1, [1,2,3],'REACCION_BASE' ], [ datos_2, [1,2,3],'DESPLAZAMIENTO_NODO2' ], [ datos_3, [1,2,3],'DESPLAZAMIENTO_NODO1' ] ]

            if serie_total:
                self.construye_figura(serie_total=serie_total, recoder='*')


        
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
                
                plt.ylabel('Desplazamiento H. (m)')
                plt.xlabel('Tiempo (s)')                       
                plt.show()
        else :
            
            from PIL import Image
            import matplotlib.image as mpimg
            import os
            from pathlib import Path
            
            BASE_DIR = Path(__file__).resolve().parent.parent

            IMAGEN_COLUMNA = os.path.join(os.path.join(BASE_DIR,'modelos'), 'columna.jpeg')

            pil_img = Image.open(IMAGEN_COLUMNA)
            

            num_subplots = 0
            
            for elemento in serie_total:
                num_subplots += len(elemento[1]) 
            
            if num_subplots % 2 != 0:
                num_subplots+=1
                        
            fig = plt.figure()


            FILAS = 3
            COLUMNAS = 4
            ax1 = plt.subplot2grid((FILAS, COLUMNAS), (0, 0), rowspan=2)
            ax1.imshow(pil_img) 
            ax1.axes.xaxis.set_visible(False)
            ax1.axes.yaxis.set_visible(False)

            ax2 = plt.subplot2grid((FILAS, COLUMNAS), (2, 0))
            ax2.axes.xaxis.set_visible(False)
            ax2.axes.yaxis.set_visible(False)
            ax2.spines['top'].set_visible(False)
            ax2.spines['right'].set_visible(False)
            ax2.spines['bottom'].set_visible(False)
            ax2.spines['left'].set_visible(False)  

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
            
                if calculo_var == 'REACCION_BASE':  
                    
                    ax3.plot(array_mek[:,0], array_mek[:,serie[0]],'tab:orange')
                    ax3.set_title(f" {calculo_var}\nGRADO DE LIBERTAD # {serie[0]}", fontsize=7)
                    ax3.set_xlabel('Tiempo (s)')
                    ax3.set_ylabel('Desplazamiento H. (m)')
                                        
                    ax4.plot(array_mek[:,0], array_mek[:,serie[1]],'tab:orange')
                    ax4.set_title(f" {calculo_var}\nGRADO DE LIBERTAD # {serie[1]}", fontsize=7)
                    ax4.set_xlabel('Tiempo (s)')
                    ax4.set_ylabel('Desplazamiento H. (m)')
                    
                    ax5.plot(array_mek[:,0], array_mek[:,serie[2]],'tab:orange')
                    ax5.set_title(f" {calculo_var}\nGRADO DE LIBERTAD # {serie[2]}", fontsize=7)
                    ax5.set_xlabel('Tiempo (s)')
                    ax5.set_ylabel('Desplazamiento H. (m)')

                elif calculo_var == 'DESPLAZAMIENTO_NODO2':
                                        
                    ax6.plot(array_mek[:,0], array_mek[:,serie[0]],'tab:red')
                    ax6.set_title(f" {calculo_var}\nGRADO DE LIBERTAD # {serie[0]}", fontsize=7)
                    ax6.set_xlabel('Tiempo (s)')
                    ax6.set_ylabel('Desplazamiento H. (m)')
                    
                    ax7.plot(array_mek[:,0], array_mek[:,serie[1]],'tab:red')
                    ax7.set_title(f" {calculo_var}\nGRADO DE LIBERTAD # {serie[1]}", fontsize=7)
                    ax7.set_xlabel('Tiempo (s)')
                    ax7.set_ylabel('Desplazamiento H. (m)')

                    ax8.plot(array_mek[:,0], array_mek[:,serie[2]],'tab:red')
                    ax8.set_title(f" {calculo_var}\nGRADO DE LIBERTAD # {serie[2]}", fontsize=7)
                    ax8.set_xlabel('Tiempo (s)')
                    ax8.set_ylabel('Desplazamiento H. (m)')      

                elif calculo_var == 'DESPLAZAMIENTO_NODO1':

                    ax9.plot(array_mek[:,0], array_mek[:,serie[0]],'tab:blue')
                    ax9.set_title(f" {calculo_var}\nGRADO DE LIBERTAD # {serie[0]}", fontsize=7)
                    ax9.set_xlabel('Tiempo (s)')
                    ax9.set_ylabel('Desplazamiento H. (m)')
                    
                    ax10.plot(array_mek[:,0], array_mek[:,serie[1]],'tab:blue')
                    ax10.set_title(f" {calculo_var}\nGRADO DE LIBERTAD # {serie[1]}", fontsize=7)
                    ax10.set_xlabel('Tiempo (s)')
                    ax10.set_ylabel('Desplazamiento H. (m)')

                    ax11.plot(array_mek[:,0], array_mek[:,serie[2]],'tab:blue')
                    ax11.set_title(f" {calculo_var}\nGRADO DE LIBERTAD # {serie[2]}", fontsize=7)
                    ax11.set_xlabel('Tiempo (s)')
                    ax11.set_ylabel('Desplazamiento H. (m)')     

                else:
                    pass         
                        
            fig.suptitle("MODELOS MATEMÁTICOS: COLUMNA AISLADOR")                             
            fig.subplots_adjust(left=0.05, bottom=0.076, right=0.924, top=0.88, wspace=0.562, hspace=0.745)
    
        mng = plt.get_current_fig_manager()
        mng.set_window_title("MODELOS MATEMÁTICOS UTPL")
        mng.resize(*mng.window.maxsize())                
        plt.show() 

    def mostrar_resultados(self, lista_mensajes=[],ax=None):
                
        fig, ax1 = plt.subplots(1, 1, constrained_layout=True, sharey=True)       
        ax1.set_xlabel('Tiempo (s)')
        ax1.set_ylabel('Desplazamiento H. (m)')   
                        
        resultado = ''
        
        for index, cadena in enumerate(lista_mensajes):           
            resultado += cadena+'\n'
                        
        ax1.axes.xaxis.set_visible(False)
        ax1.axes.yaxis.set_visible(False)
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['bottom'].set_visible(False)
        ax1.spines['left'].set_visible(False)
        ax1.text(0.038, 0.75, resultado , horizontalalignment='left', verticalalignment='center_baseline', transform=ax1.transAxes, bbox={'facecolor': 'orange', 'alpha': 0.25, 'pad': 5}, fontsize=9)        
        fig.suptitle('\n\nResultados', fontsize=16)

        plt.show()
