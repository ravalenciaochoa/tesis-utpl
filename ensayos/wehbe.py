import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt
from openseespy.postprocessing.Get_Rendering import plot_model
from PIL import Image
import matplotlib.image as mpimg
import os
from pathlib import Path


class Wehbe():

    
    FC=27200
    REC=0.028        
    SUP_INF=4
    DER_INF=5
    E_ELASTIC=200000000

    #RECODERS
    FILE_DESPLAZAMIENTO = "DESPLAZAMIENTO_1.out"  # DESPLAZAMIENTO
    FILE_REAC_BASE = "REAC_BASE.out"            # REACCION_BASE
    FILE_INTER_ENS = "INTERSECCION_ENSA_FISI.out"  # INTERSECCION DE COMPARACION DE LOS ENSAYOS REALIZADOS

    def __init__(self,nombre='BasicBuilder',numero_dimensiones=2, numero_dots=3, data_file='Wehbe et al. 1998, A1.txt', altura_columna=2.335, diametro_columna=19.1, tipo_concreto=1, tipo_refuerzo=2 ,dim_travs_colum_altura=0.61, dim_trav_colum_base=0.38, data_file_lab='Whebe.txt', calculo='*'):
                
        self.nombre = nombre
        self.numero_dimensiones = numero_dimensiones
        self.numero_dots = numero_dots                
        self.data_file = data_file
        self.tipo_concreto = tipo_concreto
        self.tipo_refuerzo = tipo_refuerzo

        #Datos de la Estructura
        self.altura_columna = altura_columna
        self.diametro_columna = diametro_columna
        self.dim_travs_colum_altura = dim_travs_colum_altura
        self.dim_trav_colum_base = dim_trav_colum_base
        self.data_file_lab = data_file_lab
        self.calculo = calculo
        
        ops.wipe()
        ops.model(self.nombre, "-ndm", self.numero_dimensiones, "-ndf", self.numero_dots)

    def construir_modelo(self):
                                                
        
        area_fiber=((np.pi*self.diametro_columna**2)/4)/1000**2

        #Definir los nodos de la estructura
        ops.node(1, 0, 0)
        ops.node(2, 0.58375, 0)
        ops.node(3, 1.1675, 0)
        ops.node(4, 1.75125, 0)
        ops.node(5, 2.335, 0)

        #Restricciones    
        ops.fix(1, 1, 1, 1)

        # Definir los materiales para la columna
        # ------------------------------------------
        # CONCRET0  
        
        if self.tipo_concreto == 1:
            # Concreto no confinado
            fpc=-self.FC
            epsc0=-0.002
            fpcu=0.8*self.FC
            epsU=-0.006
            ops.uniaxialMaterial("Concrete01", 2 , -fpc , -epsc0 , -fpcu , -epsU)

            #Concreto Confinado
            fpcc=self.FC+5521.13
            ecc=4.02982911046149E-03
            fcu=0.8*fpcc
            ecu=1.13249292446053E-02
            ops.uniaxialMaterial("Concrete01", 1, -fpcc , -ecc , -fcu , -ecu)

        elif self.tipo_concreto == 2:
            # Concreto inconfinado
            #fpc=-self.FC
            epsc0=-0.002
            fpcu=0.8*self.FC
            epsU=-0.006
            lambdac=0.5
            ft=(0.62*(self.FC)**0.5)
            ets=0.1*(2*(self.FC)/epsc0)

            ops.uniaxialMaterial("Concrete02", 2 , -self.FC , -epsc0 , -fpcu , -epsU , lambdac , ft , ets)

            # Concreto confinado

            #Columnas

            fcc=self.FC+5521.13
            epscc0=4.02982911046149E-03
            fpccu=0.8*self.FC
            epscu=1.13249292446053E-02

            ops.uniaxialMaterial("Concrete02", 1 , -fcc , -epscc0 , -fpccu , -epscu , lambdac , ft , ets)
            
        else:
            
            print("No se encontro el concreto a utilizar")

        #Material Elástico del Concreto
        
        ops.uniaxialMaterial("Elastic", 3, self.E_ELASTIC)

        #Acero de refuerzo
        
        if self.tipo_refuerzo==1:
            FY=448000
            E0=200000000
            b=1.44742225859247E-02
            a1=0.04
            a2=1
            a3=0.04
            a4=1
            ops.uniaxialMaterial("Steel01", 4, FY, E0, b, a1, a2, a3, a4)
            
        elif self.tipo_refuerzo==2:
            FY=448000
            E0=200000000
            b=1.44742225859247E-02
            R0=15
            CR1=0.925
            CR2=0.15
            a1=0.04
            a2=1
            a3=0.04
            a4=1
            ops.uniaxialMaterial("Steel02", 4, FY, E0, b, R0, CR1, CR2, a1, a2, a3, a4)

        else: 
            print("No se encontro el acero de refuerzo a utilizar")

        #Definir las etiquetas de secciones
        col_sec=1

        #Creacion de las secciones
        ops.section("Fiber", col_sec)

        y0=(self.dim_travs_colum_altura/2)-self.REC
        z0=(self.dim_trav_colum_base/2)-self.REC

        #Fibra del Concreto Confinado
        eti_quad_conf=1
        num_sub_div_ij=10
        num_sub_div_jk=10
        crds_i=[y0, z0]
        crds_j=[y0, -z0]
        crds_k=[-y0, -z0]
        crds_l=[-y0, z0]
        ops.patch("quad", eti_quad_conf, num_sub_div_ij, num_sub_div_jk, *crds_i , *crds_j , *crds_k , *crds_l)

        #Fibra del Concreto no Confinado
        y1=self.dim_travs_colum_altura/2
        z1=self.dim_trav_colum_base/2

        eti_quad_no_conf=2
        num_sub_div_ij=1
        num_sub_div_jk=4
        crds_i=[y1, z1]
        crds_j=[-y1, z1]
        crds_k=[y1, z0]
        crds_l=[-y1, z0]
        ops.patch("quad", eti_quad_no_conf, num_sub_div_ij, num_sub_div_jk, * crds_i , * crds_j , * crds_k , * crds_l)

        num_sub_div_ij=1
        num_sub_div_jk=4
        crds_i=[y1, -z1]
        crds_j=[-y1, -z1]
        crds_k=[-y1, -z0]
        crds_l=[y1, -z0]
        ops.patch("quad", eti_quad_no_conf, num_sub_div_ij, num_sub_div_jk, * crds_i , * crds_j , * crds_k , * crds_l)

        num_sub_div_ij=4
        num_sub_div_jk=1
        crds_i=[y1, z0]
        crds_j=[y1, -z0]
        crds_k=[y0, -z0]
        crds_l=[y0, z0]
        ops.patch("quad", eti_quad_no_conf, num_sub_div_ij, num_sub_div_jk, * crds_i , * crds_j , * crds_k , * crds_l)

        num_sub_div_ij=4
        num_sub_div_jk=1
        crds_i=[-y1, -z0]
        crds_j=[-y1, z0]
        crds_k=[-y0, z0]
        crds_l=[-y0, -z0]
        ops.patch("quad", eti_quad_no_conf, num_sub_div_ij, num_sub_div_jk, * crds_i , * crds_j , * crds_k , * crds_l)

        #Definir el acero longitudinal
        mat_tag=4
        #Definir acero longitudinal superior
        num_fiber=self.SUP_INF
        start=[y0, z0]
        end=[y0, -z0]
        ops.layer("straight", mat_tag, num_fiber, area_fiber, *start, *end)

        #Definir acero longitudinal inferior
        num_fiber=self.SUP_INF
        start=[-y0, z0]
        end=[-y0, -z0]
        ops.layer("straight", mat_tag, num_fiber, area_fiber, *start, *end)

        #Definir acero longitudinal derecha
        num_fiber=self.DER_INF
        start=[y0, z0]
        end=[-y0, z0]
        ops.layer("straight", mat_tag, num_fiber, area_fiber, *start, *end)

        #Definir acero longitudinal izquierda
        num_fiber=self.DER_INF
        start=[y0, -z0]
        end=[-y0, -z0]
        ops.layer("straight", mat_tag, num_fiber, area_fiber, *start, *end)

        #Definir agregador de sección 
        eti_agr=2
        ops.section("Aggregator", eti_agr, 3, "Vy", "-section", col_sec)

        #Definir la transformacion de coordenadas
        eti_transf=1
        ops.geomTransf("Linear", eti_transf)

        #Definir los elementos de la columna
        eti_beam=1
        N=10
        ops.beamIntegration("Lobatto", eti_beam, eti_agr, N )
        ele_tag=1
        ele_nodes=[1, 2]
        ops.element("dispBeamColumn", ele_tag, * ele_nodes, eti_transf, eti_beam)

        ele_tag=2
        ele_nodes=[2, 3]
        ops.element("dispBeamColumn", ele_tag, * ele_nodes, eti_transf, eti_beam)

        ele_tag=3
        ele_nodes=[3, 4]
        ops.element("dispBeamColumn", ele_tag, * ele_nodes, eti_transf, eti_beam)

        ele_tag=4
        ele_nodes=[4, 5]
        ops.element("dispBeamColumn", ele_tag, * ele_nodes, eti_transf, eti_beam)

        nodo=5
        dof=2

        #Definir el recorder
        ops.recorder("Node", "-file", self.FILE_INTER_ENS, "-time", "-node", nodo, "-dof", dof, "disp")
        ops.recorder("Node", "-file", self.FILE_DESPLAZAMIENTO, "-time", "-node", 2, "-dof", 1, 2, 3, "disp")
        ops.recorder("Node", "-file", self.FILE_REAC_BASE, "-time", "-node", 1, "-dof", 1, 2, 3, "reaction")
        
        #Definir los patrones de carga
        ops.timeSeries("Constant", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(nodo, -615, 0, 0)

        # Definir los parametros de analisis
        ops.integrator("LoadControl", 0)
        ops.system("SparseGeneral", "-piv")
        ops.test("NormDispIncr", 1.0e-5, 2000)
        ops.constraints("Plain")
        ops.numberer("Plain")
        ops.algorithm("KrylovNewton")
        ops.analysis("Static")

        ok=ops.analyze(1)
        if ok==0:
            print("CARGA AXIAL APLICADA")
        else:
            print("ERROR AL APLICAR LA CARGA AXIAL")
        ops.wipeAnalysis()

        if ok==0:
            # Definir patron de carga
            ops.timeSeries("Linear", 2)
            ops.pattern("Plain", 2, 2)
            ops.load(nodo, 0.0, 1.5, 0.0)
            
            #Definir los nodos de la estructura
            nodos=np.loadtxt(self.data_file)
            

            for i in range(len(nodos)):
                x=float(nodos[i])

                ops.integrator("DisplacementControl", nodo, dof, x)
                ops.analysis("Static")
                ops.analyze(1)


                ops.wipeAnalysis()
        ops.remove("recorders")

        #Personalizar la carga de los datos
        
        serie_total= []        

        if self.calculo == 'REACCION_BASE' :                
            array_mek=np.loadtxt(self.FILE_REAC_BASE)
            self.construye_figura(array_mek,serie=[1,2,3], recoder=self.calculo)                            


        elif self.calculo == 'DESPLAZAMIENTO_1':            
            array_mek=np.loadtxt(self.FILE_DESPLAZAMIENTO)
            self.construye_figura(array_mek,serie=[1,2,3], recoder=self.calculo)

        else:            
            datos_1=np.loadtxt(self.FILE_REAC_BASE)  
            datos_2=np.loadtxt(self.FILE_INTER_ENS)                                              
            datos_3=np.loadtxt(self.FILE_DESPLAZAMIENTO)

            serie_total = [[ datos_1, [1,2,3],'REACCION_BASE' ], [ datos_2, [1], 'INTERSECCION' ], [ datos_3, [1,2,3],'DESPLAZAMIENTO' ], ]
            #self.construye_figura(serie_total=serie_total, recoder='*')
        
        nodo=[1, 2]
        dof=[1, 2, 3]

        nodo1=nodo[0]
        nodo2=nodo[1]

        gra_liber1=dof[0]
        gra_liber2=dof[1]
        gra_liber3=dof[2]

        #Resultados

        RESULTADOS = []

        RESULTADOS.append(str("El desplazamiento en el nodo: " + str(nodo2)+ " para el grado de libertad: "+str(gra_liber1)+" es: "+str(ops.nodeDisp(2, 1))+" m") )        
        RESULTADOS.append(str("El desplazamiento en el nodo: " + str(nodo2)+ " para el grado de libertad: "+str(gra_liber2)+" es: "+str(ops.nodeDisp(2, 2))+" m") )        
        RESULTADOS.append(str("El desplazamiento en el nodo: " + str(nodo2)+ " para el grado de libertad: "+str(gra_liber3)+" es: "+str(ops.nodeDisp(2, 3))+" m") )        

        RESULTADOS.append(str("La reaccción en el nodo: " + str(nodo1)+ " para el grado de libertad: "+str(gra_liber1)+" es: "+str(ops.nodeReaction(1, 1))+" N" ))
        RESULTADOS.append(str("La reaccción en el nodo: " + str(nodo1)+ " para el grado de libertad: "+str(gra_liber2)+" es: "+str(ops.nodeReaction(1, 2))+" N") )
        RESULTADOS.append(str("La reaccción en el nodo: " + str(nodo1)+ " para el grado de libertad: "+str(gra_liber3)+" es: "+str(ops.nodeReaction(1, 3))+" N") )

        RESULTADOS.append(str("La reaccción en el elemento: "+ str(nodo1)+ " para el grado de libertad:"+ str(gra_liber1)+" es:"+str(ops.nodeReaction(1, 1))+" N"))

        if serie_total:
            self.construye_figura(serie_total=serie_total, recoder='*', resultados=RESULTADOS)

        self.mostrar_resultados(lista_mensajes=RESULTADOS) 


        #Este el el modelo realizado en el laboratorio                
        self.construye_figura([array_mek, array_mek_lab])


    def construye_figura(self, array_mek=None, serie=[], recoder=None, serie_total=[], resultados=[]):
        
        """Construye una figura basado en un array de datos, una serie de gráficos        
        """
                
        titulo = None        
        
        if recoder != '*':

            ventana_principal = plt.figure(0)
            ventana_principal.canvas.set_window_title(recoder)

            for i in serie:
                plt.subplots_adjust(hspace=0.70)
                plt.subplot(2,2,i)
                titulo = "RESPUESTA : GRADO DE LIBERTAD #{}".format(i)            
                plt.title(titulo)                        
                plt.xlabel('Tiempo (seg)')
                plt.ylabel('Desplazamiento (m)')
                print ("Normal :",i)
                plt.plot(array_mek[:,0], array_mek[:,i],label="Resultado")
                plt.legend()                           
        else :
            
            from PIL import Image
            import matplotlib.image as mpimg
            import os
            from pathlib import Path
            
            BASE_DIR = Path(__file__).resolve().parent.parent

            IMAGEN_COLUMNA = os.path.join(os.path.join(BASE_DIR,'modelos'), 'Carga_ciclica_01.jpeg')

            pil_img = Image.open(IMAGEN_COLUMNA)
            

            num_subplots = 0
            
            for elemento in serie_total:
                num_subplots += len(elemento[1]) 
            
            if num_subplots % 2 != 0:
                num_subplots+=1
                        
            fig = plt.figure()


            FILAS = 3
            COLUMNAS = 3
            ax1 = plt.subplot2grid((FILAS, COLUMNAS), (0, 0), rowspan=2)
            ax2 = plt.subplot2grid((FILAS, COLUMNAS), (2, 0))
            ax3 = plt.subplot2grid((FILAS, COLUMNAS), (0, 1))
            ax4 = plt.subplot2grid((FILAS, COLUMNAS), (1, 1))
            ax5 = plt.subplot2grid((FILAS, COLUMNAS), (2, 1))
            ax6 = plt.subplot2grid((FILAS, COLUMNAS), (0, 2))
            ax7 = plt.subplot2grid((FILAS, COLUMNAS), (1, 2))
            ax8 = plt.subplot2grid((FILAS, COLUMNAS), (2, 2))
            
            #Para controlar el subplot de Resultados se visualice
            if COLUMNAS == 4: 
                ax9 = plt.subplot2grid((3,4), (0, 3), rowspan=3)
            

            for columna, elemento in enumerate(serie_total):                
                
                array_mek = elemento[0]
                serie = elemento[1]
                calculo_var = elemento[2]
            
                if calculo_var == 'REACCION_BASE':  
                    
                    ax3.plot(array_mek[:,0], array_mek[:,serie[0]],'tab:orange')
                    ax3.set_title(f" {calculo_var}\nGRADO DE LIBERTAD # {serie[0]}", fontsize=8)
                    ax3.set_xlabel('Fuerza (kN)')
                    ax3.set_ylabel('Desplazamiento (mm)')
                                        
                    ax4.plot(array_mek[:,0], array_mek[:,serie[1]],'tab:orange')
                    ax4.set_title(f" {calculo_var}\nGRADO DE LIBERTAD # {serie[1]}", fontsize=8)
                    ax4.set_xlabel('Fuerza (kN)')
                    ax4.set_ylabel('Desplazamiento (mm)')
                    
                    ax5.plot(array_mek[:,0], array_mek[:,serie[2]],'tab:orange')
                    ax5.set_title(f" {calculo_var}\nGRADO DE LIBERTAD # {serie[2]}", fontsize=8)
                    ax5.set_xlabel('Fuerza (kN)')
                    ax5.set_ylabel('Desplazamiento (mm)')

                elif calculo_var == 'INTERSECCION':
                                    
                    ax1.imshow(pil_img) 
                    ax1.set_title(f" ENSAYO POR CARGA CÍCLICA", fontsize=8)
                    ax1.axes.xaxis.set_visible(False)
                    ax1.axes.yaxis.set_visible(False)

                    #ax4.plot(array_mek[0][:,1], array_mek[0][:,0],'tab:orange')
                    ens_opensees=np.loadtxt(self.FILE_INTER_ENS)
                    ens_lab=np.loadtxt(self.data_file_lab)
                    ax2.plot(ens_opensees[:,1], ens_opensees[:,0],'tab:green')
                    ax2.plot(ens_lab[:,0], ens_lab[:,1],'tab:orange')
                    ax2.set_title(f" INTERSECCIÓN DE ENSAYOS", fontsize=8)
                    ax2.set_xlabel('Desplazamiento (m)')
                    ax2.set_ylabel('Fuerza (kN)')


                else:

                    ax6.plot(array_mek[:,0], array_mek[:,serie[0]],'tab:red')
                    ax6.set_title(f" {calculo_var}\nGRADO DE LIBERTAD # {serie[0]}", fontsize=8)
                    ax6.set_xlabel('Fuerza (kN)')
                    ax6.set_ylabel('Desplazamiento (m)')
                    
                    ax7.plot(array_mek[:,0], array_mek[:,serie[1]],'tab:red')
                    ax7.set_title(f" {calculo_var}\nGRADO DE LIBERTAD # {serie[1]}", fontsize=8)
                    ax7.set_xlabel('Fuerza (kN)')
                    ax7.set_ylabel('Desplazamiento (m)')

                    ax8.plot(array_mek[:,0], array_mek[:,serie[2]],'tab:red')
                    ax8.set_title(f" {calculo_var}\nGRADO DE LIBERTAD # {serie[2]}", fontsize=8)
                    ax8.set_xlabel('Fuerza (kN)')
                    ax8.set_ylabel('Desplazamiento (m)')                
                        
            if COLUMNAS == 4:
                self.mostrar_resultados(lista_mensajes=resultados,ax=ax9)

            fig.suptitle("MODELOS MATEMÁTICOS: ENSAYOS FISICOS")                             
            fig.subplots_adjust(left=0.05, bottom=0.076, right=0.971, top=0.88, wspace=0.562, hspace=1)
    
        mng = plt.get_current_fig_manager()
        mng.set_window_title("MODELOS MATEMÁTICOS UTPL")
        mng.resize(*mng.window.maxsize())                
        plt.show() 

    