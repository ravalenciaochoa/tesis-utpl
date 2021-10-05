import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt

class Carga():

    FILE_DESPLAZAMIENTO_NODO2 = "DESPLAZAMIENTO_NODO2.out"          # DESPLAZAMIENTO N2
    FILE_DESPLAZAMIENTO_NODO1 = "DESPLAZAMIENTO_NODO1.out"          # DESPLAZAMIENTO N1
    FILE_REAC_BASE = "REAC_BASE.out"                                # REACCION_BASE      
    FILE_FUERZAS_COL = "FUERZAS_COL.out"                            # FUERZAS COL
    FILE_FUERZAS_COL_N1 = "FUERZAS_COL_N1.out"                      # FUERZAS COL_N1
    FILE_DEFORMACION_COL = "DEFORMACION_COL.out"                    # DEFORMACION_COL
    FILE_FORCE_COLUMNA = "Data/ForceColSecnumIntgrPts.out"          # FORCE_COLUMNA      
    FILE_DEFORMACION_COLUMNA = "Data/DefoColSec$numIntgrPts.out"    # DEFORMACION_COLUMNA

    GRAVEDAD=9.81

    def __init__(self,nombre='BasicBuilder',numero_dimensiones=2, numero_dots=3, altura_columna=12, fuerza_aplicada=15000 , resistencia_fc=21000 ,dim_travs_colum_altura=1.5, dim_trav_colum_base=1.5, calculo='*'):
        
        self.numero_dimensiones = numero_dimensiones
        self.numero_dots = numero_dots
        self.nombre = nombre
        self.dim_travs_colum_altura = dim_travs_colum_altura
        self.dim_trav_colum_base = dim_trav_colum_base
        self.calculo = calculo
        self.altura_columna = altura_columna
        self.fuerza_aplicada = fuerza_aplicada
        self.resistencia_fc = resistencia_fc

        ops.wipe()
        #Creación del modelo 
        ops.model(self.nombre, "-ndm", self.numero_dimensiones, "-ndf", self.numero_dots)

    def construir_modelo(self):

        #Datos de la estructura
                                                
        pes_esp_hor = 23.56
        peso_ele = pes_esp_hor * self.dim_travs_colum_altura * self.dim_trav_colum_base *self.altura_columna
        peso_estru = peso_ele + self.fuerza_aplicada
        
        masa = peso_estru/self.GRAVEDAD

        area = self.dim_travs_colum_altura*self.dim_trav_colum_base
        
        E=20000000

        #Definir los nodos de la estructura
        ops.node(1, 0, 0)
        ops.node(2, 0, self.altura_columna)

        #Definir las restricciones
        ops.fix(1, 1, 1, 1)

        #Definir la masa de la estructura
        ops.mass(2, 0, 330000, 0)
        
        ea_col=E*area
        my_col=130000
        phy_col=0.65e-4
        el_col_crack=my_col/phy_col
        b=0.01

        ops.uniaxialMaterial("Steel01", 2, my_col, el_col_crack, b)
        ops.uniaxialMaterial("Elastic", 3, ea_col)
        ops.section("Aggregator", 1, 3, "P", 2, "Mz")

        #Definir la transformacion de coordenadas
        ops.geomTransf("Linear", 1)
        ops.element("nonlinearBeamColumn", 1, 1, 2, 5, 1, 1)

        #Definir los recorders

        ops.recorder("Node", "-file", self.FILE_DESPLAZAMIENTO_NODO2, "-time", "-node", 2, "-dof", 1, 2, 3, "disp")
        ops.recorder("Node", "-file", self.FILE_DESPLAZAMIENTO_NODO1, "-time", "-node", 1, "-dof", 1, 2, 3, "disp")
        ops.recorder("Node", "-file", self.FILE_REAC_BASE, "-time", "-node", 1, "-dof", 1, 2, 3, "reaction")
        
        ops.recorder("Element", "-file", self.FILE_FUERZAS_COL, "-time", "-ele", 2, "globalForce")
        ops.recorder("Element", "-file", self.FILE_FUERZAS_COL_N1, "-time", "-ele", 1, "section", 1, "force")
        ops.recorder("Element", "-file", self.FILE_DEFORMACION_COL, "-time", "-ele", 1, "section", 1, "deformation")
        ops.recorder("Element", "-file", self.FILE_FORCE_COLUMNA, "-time", "-ele", 1, "section", 5, "force")
        ops.recorder("Element", "-file", self.FILE_DEFORMACION_COLUMNA, "-time", "-ele", 1, "section", 1, "deformation")

        #Definir las actuantes sobre la columna
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(2, 0, -peso_estru, 0)

        #Parametros de analisis estatico
        ops.constraints("Plain")
        ops.numberer("Plain")
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

        if ok == 0:
            # Definir patron de carga
            ops.timeSeries("Linear", 2)
            ops.pattern("Plain", 2, 2)
            ops.load(2, peso_estru, 0.0, 0.0, 0.0, 0.0, 0.0)

            #Parametros de analisis estatico
            ops.constraints("Plain")
            ops.numberer("Plain")
            ops.system("BandGeneral")
            ops.test("EnergyIncr", 1.0e-8, 10, 0)
            ops.algorithm("Newton")
            max_disp=0.05*self.altura_columna
            max_steps=10

            delta_disp=max_disp/max_steps
            ops.integrator("DisplacementControl", 2, 1, delta_disp)
            ops.analysis("Static")
            ops.analyze((0.01*self.altura_columna)/(0.001*self.altura_columna))

            while ops.nodeDisp(2, 1) <= max_disp:

                for n in [1, 2, 4, 10]:
                    ops.integrator("DisplacementControl", 2, 1, delta_disp/n)
                    ok=ops.analyze(n)

                    if ok !=0:
                        ops.test("NormDispIncr", 1.0e-5, 5000)
                        ops.algorithm("Newton", "-initial")
                        ok=ops.analyze(n)
                        ops.test("EnergyIncr", 1.0e-6, 100)
                    if ok !=0:
                        ops.test("NormDispIncr", 1.0e-4, 10000)
                        ops.algorithm("Newton", "-initial")
                        ok=ops.analyze(n)
                        ops.test("EnergyIncr", 1.0e-6, 100)
                    if ok !=0:
                        ops.test("NormDispIncr", 1.0e-3, 10000)
                        ops.algorithm("Newton", "-initial")
                        ok=ops.analyze(n)
                        ops.test("EnergyIncr", 1.0e-6, 100)
                    if ok !=0:
                        ops.test("NormDispIncr", 1.0e-2, 10000)
                        ops.algorithm("Newton", "-initial")                        
                        ok=ops.analyze(n)
                        ops.test("EnergyIncr", 1.0e-6, 100)
                    
                    if ok ==0:
                        break

                if ok != 0:
                    print("EL ANÁLISIS DE DESPLAZAMIENTO FALLÓ: ", ops.nodeDisp(2, 2))
                    break

            if ok == 0:
                print("EL ANÁLISIS DE DESPLAZAMIENTO TERMINÓ EXITOSAMENTE")
                
            ops.wipeAnalysis()
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
        
        elif self.calculo == 'FUERZAS_COL':            
            array_mek=np.loadtxt(self.FILE_FUERZAS_COL)           
            self.construye_figura(array_mek,serie=[1], recoder=self.calculo)

        elif self.calculo == 'FUERZAS_COL_N1':            
            array_mek=np.loadtxt(self.FILE_FUERZAS_COL_N1)
            self.construye_figura(array_mek,serie=[1], recoder=self.calculo)
        
        elif self.calculo == 'DEFORMACION_COL':            
            array_mek=np.loadtxt(self.FILE_DEFORMACION_COL)
            self.construye_figura(array_mek,serie=[1], recoder=self.calculo)

        elif self.calculo == 'FORCE_COLUMNA':            
            array_mek=np.loadtxt(self.FILE_FORCE_COLUMNA)
            self.construye_figura(array_mek,serie=[1], recoder=self.calculo)

        elif self.calculo == 'DEFORMACION_COLUMNA':            
            array_mek=np.loadtxt(self.FILE_DEFORMACION_COLUMNA)
            self.construye_figura(array_mek,serie=[1], recoder=self.calculo)

        else:            
            datos_1=np.loadtxt(self.FILE_REAC_BASE)
            #self.construye_figura(array_mek,serie=[1,2,3], recoder='REACCION_BASE')

            datos_2=np.loadtxt(self.FILE_DESPLAZAMIENTO_NODO2)
            #self.construye_figura(array_mek,serie=[1,2,3], recoder='DESPLAZAMIENTO_NODO2')

            datos_3=np.loadtxt(self.FILE_DESPLAZAMIENTO_NODO1)
            #self.construye_figura(array_mek,serie=[1,2,3], recoder='DESPLAZAMIENTO_NODO1')


            """array_mek=np.loadtxt(self.FILE_FUERZAS_COL)
            #self.construye_figura(array_mek,serie=[1], recoder='FUERZAS_COL')

            array_mek=np.loadtxt(self.FILE_FUERZAS_COL_N1)
            #self.construye_figura(array_mek,serie=[1], recoder='FUERZAS_COL_N1')

            array_mek=np.loadtxt(self.FILE_DEFORMACION_COL)
            #self.construye_figura(array_mek,serie=[1], recoder='DEFORMACION_COL')

            array_mek=np.loadtxt(self.FILE_FORCE_COLUMNA)
            #self.construye_figura(array_mek,serie=[1], recoder='FORCE_COLUMNA')

            array_mek=np.loadtxt(self.FILE_DEFORMACION_COLUMNA)
            #self.construye_figura(array_mek,serie=[1], recoder='DEFORMACION_COLUMNA')"""

            serie_total = [ [ datos_1, [1,2,3],'REACCION_BASE' ], [ datos_2, [1,2,3],'DESPLAZAMIENTO_NODO2' ], [ datos_3, [1,2,3],'DESPLAZAMIENTO_NODO1' ] ]
        
        Num_modos=1
        OmegaSq=ops.eigen("-fullGenLapack", Num_modos)
        OmegaSq=np.array(OmegaSq)
        Omega=OmegaSq**0.5
        Periodo=(2*np.pi)/Omega

        for i in range(Num_modos):
            print("Modo", i+1, "T= ", Periodo[i], "s")

        nodo=[1, 2]
        dof=[1, 2, 3]

        nodo1=nodo[0]
        nodo2=nodo[1]

        gra_liber1=dof[0]
        gra_liber2=dof[1]
        gra_liber3=dof[2]

        #Resultados

        RESULTADOS = []

        RESULTADOS.append(str("El desplazamiento en el nodo " + str(nodo2) +" para el grado de libertad " + str(gra_liber1) + " es: " + str(ops.nodeDisp(2, 1))+" m" ))
        RESULTADOS.append(str("El desplazamiento en el nodo " + str(nodo2) +" para el grado de libertad " + str(gra_liber2) + " es: " + str(ops.nodeDisp(2, 2))+" m" ))
        RESULTADOS.append(str("El desplazamiento en el nodo " + str(nodo2) +" para el grado de libertad " + str(gra_liber3) + " es: " + str(ops.nodeDisp(2, 3))+" rad/s" ))

        RESULTADOS.append(str("La reaccción en el nodo "+ str(nodo1) +" para el grado de libertad " + str(gra_liber1) + " es: " + str(ops.nodeReaction(1, 1))+" kN" ))
        RESULTADOS.append(str("La reaccción en el nodo "+ str(nodo1) +" para el grado de libertad " + str(gra_liber2) + " es: " + str(ops.nodeReaction(1, 2))+" kN" ))
        RESULTADOS.append(str("La reaccción en el nodo "+ str(nodo1) +" para el grado de libertad " + str(gra_liber3) + " es: " + str(ops.nodeReaction(1, 3))+" rad/s" ))

        
        if serie_total:
            self.construye_figura(serie_total=serie_total, recoder='*', resultados=RESULTADOS)

        self.mostrar_resultados(lista_mensajes=RESULTADOS)

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
                
                plt.plot(array_mek[:,0], array_mek[:,i], '-x')
                plt.ylabel('Carga aplicada (kN)')
                plt.xlabel('Desplazamiento horizontal (mm)')                
                plt.show()
        else :
            
            from PIL import Image
            import matplotlib.image as mpimg
            import os
            from pathlib import Path
            
            BASE_DIR = Path(__file__).resolve().parent.parent

            IMAGEN_COLUMNA = os.path.join(os.path.join(BASE_DIR,'modelos'), 'Columna_con_carga_puntual.jpeg')

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
                    ax3.set_title(f" {calculo_var}\nGRADO DE LIBERTAD # {serie[0]}", fontsize=8)
                    ax3.set_xlabel('Desplazamiento Horizontal (mm)')
                    ax3.set_ylabel('Carga Aplicada  (Kn)')
                                        
                    ax4.plot(array_mek[:,0], array_mek[:,serie[1]],'tab:orange')
                    ax4.set_title(f" {calculo_var}\nGRADO DE LIBERTAD # {serie[1]}", fontsize=8)
                    ax4.set_xlabel('Desplazamiento Horizontal (mm)')
                    ax4.set_ylabel('Carga Aplicada  (Kn)')
                    
                    ax5.plot(array_mek[:,0], array_mek[:,serie[2]],'tab:orange')
                    ax5.set_title(f" {calculo_var}\nGRADO DE LIBERTAD # {serie[2]}", fontsize=8)
                    ax5.set_xlabel('Desplazamiento Horizontal (mm)')
                    ax5.set_ylabel('Carga Aplicada  (Kn)')

                elif calculo_var == 'DESPLAZAMIENTO_NODO2':
                                        
                    ax6.plot(array_mek[:,0], array_mek[:,serie[0]],'tab:red')
                    ax6.set_title(f" {calculo_var}\nGRADO DE LIBERTAD # {serie[0]}", fontsize=8)
                    ax6.set_xlabel('Desplazamiento Horizontal (mm)')
                    ax6.set_ylabel('Carga Aplicada  (Kn)')
                    
                    ax7.plot(array_mek[:,0], array_mek[:,serie[1]],'tab:red')
                    ax7.set_title(f" {calculo_var}\nGRADO DE LIBERTAD # {serie[1]}", fontsize=8)
                    ax7.set_xlabel('Desplazamiento Horizontal (mm)')
                    ax7.set_ylabel('Carga Aplicada  (Kn)')

                    ax8.plot(array_mek[:,0], array_mek[:,serie[2]],'tab:red')
                    ax8.set_title(f" {calculo_var}\nGRADO DE LIBERTAD # {serie[2]}", fontsize=8)
                    ax8.set_xlabel('Desplazamiento Horizontal (mm)')
                    ax8.set_ylabel('Carga Aplicada  (Kn)')      

                elif calculo_var == 'DESPLAZAMIENTO_NODO1':

                    ax9.plot(array_mek[:,0], array_mek[:,serie[0]],'tab:blue')
                    ax9.set_title(f" {calculo_var}\nGRADO DE LIBERTAD # {serie[0]}", fontsize=8)
                    ax9.set_xlabel('Desplazamiento Horizontal (mm)')
                    ax9.set_ylabel('Carga Aplicada  (Kn)')
                    
                    ax10.plot(array_mek[:,0], array_mek[:,serie[1]],'tab:blue')
                    ax10.set_title(f" {calculo_var}\nGRADO DE LIBERTAD # {serie[1]}", fontsize=8)
                    ax10.set_xlabel('Desplazamiento Horizontal (mm)')
                    ax10.set_ylabel('Carga Aplicada  (Kn)')

                    ax11.plot(array_mek[:,0], array_mek[:,serie[2]],'tab:blue')
                    ax11.set_title(f" {calculo_var}\nGRADO DE LIBERTAD # {serie[2]}", fontsize=8)
                    ax11.set_xlabel('Desplazamiento Horizontal (mm)')
                    ax11.set_ylabel('Carga Aplicada  (Kn)')     

                else:
                    pass         
                        
            fig.suptitle("MODELOS MATEMÁTICOS: COLUMNA DE HORMIGÓN ARMADO CON CARGA PUNTUAL APLICADA EN SENTIDO VERTICAL")                             
            fig.subplots_adjust(left=0.05, bottom=0.076, right=0.971, top=0.88, wspace=0.562, hspace=1)
    
        mng = plt.get_current_fig_manager()
        mng.set_window_title("MODELOS MATEMÁTICOS UTPL")
        mng.resize(*mng.window.maxsize())                
        plt.show() 

    def mostrar_resultados(self, lista_mensajes=[],ax=None):
                
        fig, ax1 = plt.subplots(1, 1, constrained_layout=True, sharey=True)       
        ax1.set_xlabel('Desplazamiento Horizontal (mm)')
        ax1.set_ylabel('Carga Aplicada  (Kn)')   
                        
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



        
        
        
        