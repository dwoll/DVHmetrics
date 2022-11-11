using System;
using System.Linq;
using System.Text;
using System.Windows;
using System.Collections.Generic;
using VMS.TPS.Common.Model.API;
using VMS.TPS.Common.Model.Types;
using System.IO;
using System.Windows.Media.Media3D;
using System.Windows.Media;

namespace VMS.TPS
{
  public class Script
  {
    public Script()
    {
    }

    //---------------------------------------------------------------------------------------------
    public void Execute(ScriptContext context /*, System.Windows.Window window*/)
    {
      Patient patient = context.Patient;
      if (context.StructureSet == null)
      {
        MessageBox.Show("Please, open an image with structures");
        return;
      }

      string folder = @"c:\temp\Export3D\" + MakeFilenameValid(patient.Id);
      if (!Directory.Exists(folder))
      {
        Directory.CreateDirectory(folder);
      }

      ExportStructures(context.StructureSet, folder);
      if (context.PlanSetup != null && context.PlanSetup.Dose != null)
      {
        context.PlanSetup.DoseValuePresentation = DoseValuePresentation.Absolute;
        ExportDose(context.PlanSetup.Dose, folder);
      }
      MessageBox.Show("Exported 3D data to " + folder);
    }

    //---------------------------------------------------------------------------------------------
    void ExportStructures(StructureSet ss, string folder)
    {
      foreach (Structure structure in ss.Structures)
      {
        if (!structure.HasSegment)
          continue;
        string id = structure.Id;
        string filename = MakeFilenameValid(id);
        SaveTriangleMeshToPlyFile(structure.MeshGeometry, structure.Color, folder + "\\" + filename + ".ply");
      }
    }

    //---------------------------------------------------------------------------------------------
    void ExportDose(Dose dose, string folder)
    {
      string filename = folder + "\\dose.vtk";
      SaveDoseOrImageToVTKStructurePoints(dose, null, filename);

      foreach (Isodose isodose in dose.Isodoses)
      {
        string id = isodose.Level.ToString();
        filename = MakeFilenameValid(id.Replace(" ", "_").Replace("%", "p"));
        SaveTriangleMeshToPlyFile(isodose.MeshGeometry, isodose.Color, folder + "\\" + filename + ".ply");
      }
    }

    //---------------------------------------------------------------------------------------------
    /// <summary>
    /// This method saves the given MeshGeometry3D to the given file in the PLY format
    /// also known as Polygon File Format or Stanfor Triangle Formatn
    /// </summary>
    /// <param name="mesh">Trianglemesh to export</param>
    /// <param name="col">Color of the mesh</param>
    /// <param name="outputFileName">Name of the file to write.</param>
    //---------------------------------------------------------------------------------------------
    void SaveTriangleMeshToPlyFile(MeshGeometry3D mesh, Color color, string outputFileName)
    {
      if (mesh == null)
        return;

      if (File.Exists(outputFileName))
      {
        File.SetAttributes(outputFileName, FileAttributes.Normal);
        File.Delete(outputFileName);
      }

      Point3DCollection vertexes = mesh.Positions;
      Int32Collection indexes = mesh.TriangleIndices;

      byte alpha = (byte)(0.6*255);

      using (TextWriter writer = new StreamWriter(outputFileName))
      {
        writer.WriteLine("ply");
        writer.WriteLine("format ascii 1.0");
        writer.WriteLine("element vertex " + vertexes.Count);

        writer.WriteLine("property float x");
        writer.WriteLine("property float y");
        writer.WriteLine("property float z");

        writer.WriteLine("property uchar red");
        writer.WriteLine("property uchar green");
        writer.WriteLine("property uchar blue");
        writer.WriteLine("property uchar alpha");

        writer.WriteLine("element face " + indexes.Count / 3);

        writer.WriteLine("property list uchar int vertex_indices");

        writer.WriteLine("end_header");

        foreach (Point3D v in vertexes)
        {
          writer.Write(v.X.ToString("e") + " ");
          writer.Write(v.Y.ToString("e") + " ");
          writer.Write(v.Z.ToString("e") + " ");

          writer.Write(color.R + " ");
          writer.Write(color.G + " ");
          writer.Write(color.B + " ");
          writer.Write(alpha);
          writer.WriteLine();
        }

        int i = 0;
        while (i < indexes.Count)
        {
          writer.Write("3 ");
          writer.Write(indexes[i++] + " ");
          writer.Write(indexes[i++] + " ");
          writer.Write(indexes[i++] + " ");
          writer.WriteLine();
        }
      }
    }

    //---------------------------------------------------------------------------------------------
    /// <summary>
    /// This method saves the given dose or image voxels to a file in the vtk structure points format.
    /// Note that either Dose or Image is given as parameter and the other must be null.
    /// </summary>
    /// <param name="dose">Dose to output or null</param>
    /// <param name="image">Image to output or null</param>
    /// <param name="outputFileName">>Name of the file to write</param>
    //---------------------------------------------------------------------------------------------
    public static void SaveDoseOrImageToVTKStructurePoints(Dose dose, Image image, string outputFileName)
    {
      if (File.Exists(outputFileName))
      {
        File.SetAttributes(outputFileName, FileAttributes.Normal);
        File.Delete(outputFileName);
      }

      int W, H, D;
      double sx, sy, sz;
      VVector origin, rowDirection, columnDirection;
      if (dose != null)
      {
        W = dose.XSize;
        H = dose.YSize;
        D = dose.ZSize;
        sx = dose.XRes;
        sy = dose.YRes;
        sz = dose.ZRes;
        origin = dose.Origin;
        rowDirection = dose.XDirection;
        columnDirection = dose.YDirection;
      }
      else
      {
        W = image.XSize;
        H = image.YSize;
        D = image.ZSize;
        sx = image.XRes;
        sy = image.YRes;
        sz = image.ZRes;
        origin = image.Origin;
        rowDirection = image.XDirection;
        columnDirection = image.YDirection;
      }

      using (TextWriter writer = new StreamWriter(outputFileName))
      {
        writer.WriteLine("# vtk DataFile Version 3.0");
        writer.WriteLine("vtk output");
        writer.WriteLine("ASCII");
        writer.WriteLine("DATASET STRUCTURED_POINTS");
        writer.WriteLine("DIMENSIONS " + W + " " + H + " " + D);

        int[,] buffer = new int[W, H];

        double xsign = rowDirection.x > 0 ? 1.0 : -1.0;
        double ysign = columnDirection.y > 0 ? 1.0 : -1.0;
        double zsign = GetZDirection(rowDirection, columnDirection).z > 0 ? 1.0 : -1.0;

        writer.WriteLine("ORIGIN " + origin.x.ToString() + " " + origin.y.ToString() + " " + origin.z.ToString());
        writer.WriteLine("SPACING " + sx * xsign + " " + sy * ysign + " " + sz * zsign);
        writer.WriteLine("POINT_DATA " + W * H * D);
        writer.WriteLine("SCALARS image_data unsigned_short 1");
        writer.WriteLine("LOOKUP_TABLE default");

        int maxValueForScaling = dose != null ? FindMaxValue(dose) : 0;

        for (int z = 0; z < D; z++)
        {
          if (dose != null) dose.GetVoxels(z, buffer);
          else image.GetVoxels(z, buffer);
          for (int y = 0; y < H; y++)
          {
            for (int x = 0; x < W; x++)
            {
              int value = buffer[x, y];
              UInt16 curvalue = 0;
              if (image != null)
                curvalue = (UInt16)value;
              else
                curvalue = (UInt16)((double)value * 100 / (double)maxValueForScaling);
              writer.Write(curvalue + " ");
            }
            writer.WriteLine();
          }
        }
      }
    }

    //---------------------------------------------------------------------------------------------
    /// <summary>
    /// This method finds the maximum vocel value of the dose matrix
    /// </summary>
    /// <param name="dose">Dose</param>
    /// <returns>Maximum value</returns>
    //---------------------------------------------------------------------------------------------
    static int FindMaxValue(Dose dose)
    {
      int maxValue = 0;
      int[,] buffer = new int[dose.XSize, dose.YSize];
      if (dose != null)
      {
        for (int z = 0; z < dose.ZSize; z++)
        {
          dose.GetVoxels(z, buffer);
          for (int y = 0; y < dose.YSize; y++)
          {
            for (int x = 0; x < dose.XSize; x++)
            {
              int value = buffer[x, y];
              if (value > maxValue)
                maxValue = value;
            }
          }
        }
      }
      return maxValue;
    }

    private static VVector GetZDirection(VVector a, VVector b)
    {
      // return cross product
      return new VVector(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
    }


    string MakeFilenameValid(string s)
    {
      char[] invalidChars = System.IO.Path.GetInvalidFileNameChars();
      foreach (char ch in invalidChars)
      {
        s = s.Replace(ch, '_');
      }
      return s;
    }
  }
}
