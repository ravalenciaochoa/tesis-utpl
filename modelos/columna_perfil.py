import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt
from openseespy.postprocessing.Get_Rendering import plot_model


class ColumnaPerfil():

    FILE_DESPLAZAMIENTO_NODO2 = "DESPLAZAMIENTO_NODO2.out"          # DESPLAZAMIENTO N2
    FILE_DESPLAZAMIENTO_NODO1 = "DESPLAZAMIENTO_NODO1.out"          # DESPLAZAMIENTO N1
    FILE_REAC_BASE = "REAC_BASE.out"                                # REACCION_BASE      
    FILE_FUERZAS_COL = "FUERZAS_COL.out"                            # FUERZAS COL
    FILE_FUERZAS_COL_N1 = "FUERZAS_COL_N1.out"                      # FUERZAS COL_N1
    FILE_DEFORMACION_COL = "DEFORMACION_COL.out"                    # DEFORMACION_COL
    FILE_FORCE_COLUMNA = "Data/ForceColSecnumIntgrPts.out"          # FORCE_COLUMNA      
    FILE_DEFORMACION_COLUMNA = "Data/DefoColSec$numIntgrPts.out"    # DEFORMACION_COLUMNA

    GRAVEDAD=9.81


    def __init__(self,nombre='BasicBuilder',numero_dimensiones=2, numero_dots=3, altura_columna=12, fuerza_aplicada=1036.21 , resistencia_fc=21000 ,dim_travs_colum_altura=1.5, dim_trav_colum_base=1.5, calculo='*'):
        
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
        
        peso=self.fuerza_aplicada/self.GRAVEDAD

        #Definir los nodos de la estructura
        ops.node(1, 0, 0)
        ops.node(2, 0, self.altura_columna)

        #Restricciones    
        ops.fixY(0, 1, 1, 1)

        #Definir los materiales
        fy=420000
        E=200000000
        b=0.007
        R0=20
        cr1=0.925
        cr2=0.15
        a1=0.04
        a2=1
        a3=0.04
        a4=1

        params=[R0, cr1, cr2]
        ops.uniaxialMaterial("Steel02", 1, fy, E, b, *params, a1, a2, a3, a4)

        #Definir las etiquetas de secciones
        col_sec=1

        #Creacion de las secciones
        rec=0.02
        y0=self.dim_trav_colum_base/2
        y1=y0-rec
        z0=self.dim_travs_colum_altura/2
        z1=z0-rec
        num_div_z=16
        num_div_y=16
        num_div_rec=2
        ops.section("Fiber", col_sec)

        # Definir los parches para recubriento lateral

        crdsi=[-y0, -z0]
        crdsj=[-y1, -z0]
        crdsk=[-y1, +z0]
        crdsl=[-y0, +z0]

        ops.patch("quad", 1, num_div_rec, (num_div_z+2*num_div_rec), *crdsi, *crdsj, *crdsk, *crdsl)

        crdsi=[y1, -z0]
        crdsj=[y0, -z0]
        crdsk=[y0, +z0]
        crdsl=[y1, +z0]

        ops.patch("quad", 1, num_div_rec, (num_div_z+2*num_div_rec), *crdsi, *crdsj, *crdsk, *crdsl)

        # Definir los parches para recubriento inferior y superior
        crdsi=[-y1, -z0]
        crdsj=[y1, -z0]
        crdsk=[y1, -z1]
        crdsl=[-y1, -z1]

        ops.patch("quad", 1, num_div_y, num_div_rec, *crdsi, *crdsj, *crdsk, *crdsl)

        crdsi=[-y1, -z1]
        crdsj=[y1, -z1]
        crdsk=[y1, +z0]
        crdsl=[-y1, +z0]

        ops.patch("quad", 1, num_div_y, num_div_rec, *crdsi, *crdsj, *crdsk, *crdsl)

        num_subdiv_circ=5
        num_subdiv_rad=5
        centro=[-self.dim_travs_colum_altura/2, -self.dim_trav_colum_base/2]
        rad=[1.58, 3.16]
        ang=[180, 270]
        ops.patch( "circ" , 1, num_subdiv_circ, num_subdiv_rad, * centro, * rad, * ang)

        num_subdiv_circ=5
        num_subdiv_rad=5
        centro=[self.dim_travs_colum_altura/2, -self.dim_trav_colum_base/2]
        rad=[1.58, 3.16]
        ang=[270, 360]
        ops.patch( "circ" , 1, num_subdiv_circ, num_subdiv_rad, * centro, * rad, * ang)

        num_subdiv_circ=5
        num_subdiv_rad=5
        centro=[self.dim_travs_colum_altura/2, self.dim_trav_colum_base/2]
        rad=[1.58, 3.16]
        ang=[0, 90]
        ops.patch( "circ" , 1, num_subdiv_circ, num_subdiv_rad, * centro, * rad, * ang)

        num_subdiv_circ=5
        num_subdiv_rad=5
        centro=[-self.dim_travs_colum_altura/2, self.dim_trav_colum_base/2]
        rad=[1.58, 3.16]
        ang=[90, 180]
        ops.patch( "circ" , 1, num_subdiv_circ, num_subdiv_rad, * centro, * rad, * ang)

        #Definir las etiquetas de transformacion
        col_trans=1

        #Definir los elementos de la columna
        ops.geomTransf("Corotational", col_trans)
        col_intr=1
        N=5
        ops.beamIntegration( 'Lobatto', col_intr, col_sec , N )
        ops.element("forceBeamColumn", 1, 1, 2, col_intr, col_sec, col_trans)

        #Definir los patrones de carga
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(2, 0, -self.fuerza_aplicada, 0)

        #Definir los recorders
        ops.recorder("Node", "-file", self.FILE_DESPLAZAMIENTO_NODO2, "-time", "-node", 2, "-dof", 1, 2, 3, "disp")
        ops.recorder("Node", "-file", self.FILE_DESPLAZAMIENTO_NODO1, "-time", "-node", 1, "-dof", 1, 2, 3, "disp")
        ops.recorder("Node", "-file", self.FILE_REAC_BASE, "-time", "-node", 1, "-dof", 1, 2, 3, "reaction")

        #Crear el analisis
        ops.constraints("Plain")
        ops.numberer("RCM")
        ops.system("UmfPack")
        ops.test("NormDispIncr", 1.0e-5, 10, 0)
        ops.algorithm("Newton")
        num_pasos=10
        ops.integrator("LoadControl", 1/num_pasos)
        ops.analysis("Static")
        ok=ops.analyze(num_pasos)

        if ok==0:
            print("CARGA AXIAL APLICADA")
        else:
            print("ERROR AL APLICAR LA CARGA AXIAL")
            
        ops.loadConst("-time", 0.0)

        if ok==0:
            # Patron de carga axial
            ops.timeSeries("Linear", 2)
            ops.pattern("Plain", 2, 2)
        
            ops.load(2, 0.0, 0.0, 0.0)
                
            #Parametros de analisis estatico
            id_node=2
            d_max_push=0.5
            dx_push=100
            nsteps=d_max_push/dx_push

            ops.constraints("Plain")
            ops.numberer("RCM")
            ops.system("UmfPack")
            ops.test("NormDispIncr", 1.0e-3, 500, 0)
            ops.algorithm("Newton")
            ops.integrator("DisplacementControl", id_node, 1, dx_push)
            ops.analysis("Static")
            
            while ops.nodeDisp(2, 1) <= d_max_push:
                for n in [1, 2, 4, 10]:
                    ops.integrator("DisplacementControl", id_node, 1, nsteps/n)
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
                if ok !=0:
                    print("EL ANÁLISIS DE DESPLAZAMIENTO FALLÓ: ", ops.nodeDisp(2, 1))
                    break
            if ok==0:
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




        
    def construye_figura(self, array_mek=None, serie=[], recoder=None, serie_total=[], RESULTADOS=[]):

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
                
                plt.plot(array_mek[:,0], array_mek[:,i], "-.")
                
                plt.ylabel('Carga aplicada (N)')
                plt.xlabel('Desplazamiento horizontal (m)')
                plt.show()
        else :
            
            from PIL import Image
            import matplotlib.image as mpimg
            import os
            from pathlib import Path
            
            BASE_DIR = Path(__file__).resolve().parent.parent

            IMAGEN_COLUMNA = os.path.join(os.path.join(BASE_DIR,'modelos'), 'Columna_con_perfil_metalico.jpeg')

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
                    ax3.set_xlabel('Desplazamiento Horizontal (m)')
                    ax3.set_ylabel('Carga Aplicada  (N)')
                                        
                    ax4.plot(array_mek[:,0], array_mek[:,serie[1]],'tab:orange')
                    ax4.set_title(f" {calculo_var}\nGRADO DE LIBERTAD # {serie[1]}", fontsize=8)
                    ax4.set_xlabel('Desplazamiento Horizontal (mm)')
                    ax4.set_ylabel('Carga Aplicada  (N)')
                    
                    ax5.plot(array_mek[:,0], array_mek[:,serie[2]],'tab:orange')
                    ax5.set_title(f" {calculo_var}\nGRADO DE LIBERTAD # {serie[2]}", fontsize=8)
                    ax5.set_xlabel('Desplazamiento Horizontal (mm)')
                    ax5.set_ylabel('Carga Aplicada  (N)')

                elif calculo_var == 'DESPLAZAMIENTO_NODO2':
                                        
                    ax6.plot(array_mek[:,0], array_mek[:,serie[0]],'tab:red')
                    ax6.set_title(f" {calculo_var}\nGRADO DE LIBERTAD # {serie[0]}", fontsize=8)
                    ax6.set_xlabel('Desplazamiento Horizontal (mm)')
                    ax6.set_ylabel('Carga Aplicada  (N)')
                    
                    ax7.plot(array_mek[:,0], array_mek[:,serie[1]],'tab:red')
                    ax7.set_title(f" {calculo_var}\nGRADO DE LIBERTAD # {serie[1]}", fontsize=8)
                    ax7.set_xlabel('Desplazamiento Horizontal (mm)')
                    ax7.set_ylabel('Carga Aplicada  (N)')

                    ax8.plot(array_mek[:,0], array_mek[:,serie[2]],'tab:red')
                    ax8.set_title(f" {calculo_var}\nGRADO DE LIBERTAD # {serie[2]}", fontsize=8)
                    ax8.set_xlabel('Desplazamiento Horizontal (mm)')
                    ax8.set_ylabel('Carga Aplicada  (N)')      

                elif calculo_var == 'DESPLAZAMIENTO_NODO1':

                    ax9.plot(array_mek[:,0], array_mek[:,serie[0]],'tab:blue')
                    ax9.set_title(f" {calculo_var}\nGRADO DE LIBERTAD # {serie[0]}", fontsize=8)
                    ax9.set_xlabel('Desplazamiento Horizontal (mm)')
                    ax9.set_ylabel('Carga Aplicada  (N)')
                    
                    ax10.plot(array_mek[:,0], array_mek[:,serie[1]],'tab:blue')
                    ax10.set_title(f" {calculo_var}\nGRADO DE LIBERTAD # {serie[1]}", fontsize=8)
                    ax10.set_xlabel('Desplazamiento Horizontal (mm)')
                    ax10.set_ylabel('Carga Aplicada  (N)')

                    ax11.plot(array_mek[:,0], array_mek[:,serie[2]],'tab:blue')
                    ax11.set_title(f" {calculo_var}\nGRADO DE LIBERTAD # {serie[2]}", fontsize=8)
                    ax11.set_xlabel('Desplazamiento Horizontal (mm)')
                    ax11.set_ylabel('Carga Aplicada  (N)')     

                else:
                    pass         
                        
            fig.suptitle("MODELOS MATEMÁTICOS: COLUMNA MIXTA HORMIGÓN-ACERO CON CARGA PUNTUAL")                             
            fig.subplots_adjust(left=0.05, bottom=0.121, right=0.902, top=0.88, wspace=0.562, hspace=1)
    
        mng = plt.get_current_fig_manager()
        mng.set_window_title("MODELOS MATEMÁTICOS UTPL")
        mng.resize(*mng.window.maxsize())                
        plt.show() 

    def mostrar_resultados(self, lista_mensajes=[],ax=None):
                
        fig, ax1 = plt.subplots(1, 1, constrained_layout=True, sharey=True)       
        ax1.set_xlabel('Desplazamiento Horizontal (mm)')
        ax1.set_ylabel('Carga Aplicada  (N)')   
                        
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

    