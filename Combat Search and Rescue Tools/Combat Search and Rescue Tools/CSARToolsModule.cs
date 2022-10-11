using ArcGIS.Desktop.Core;
using ArcGIS.Desktop.Core.Events;
using ArcGIS.Desktop.Framework;
using ArcGIS.Desktop.Framework.Contracts;
using System;
using System.IO;

namespace Combat_Search_and_Rescue_Tools
{
    internal class CSARToolsModule : Module
    {
        private static CSARToolsModule _this = null;
        internal const string AddInID = "{771ace14-71e9-4a6f-a17d-7aaf70afcc59}";
        private static readonly string _installPath = Path.Combine(Environment.GetFolderPath(Environment.SpecialFolder.LocalApplicationData), string.Format(@"ESRI\ArcGISPro\Toolboxes\{0}\toolboxes", AddInID));


        /// <summary>
        /// Retrieve the singleton instance to this module here
        /// </summary>
        public static CSARToolsModule Current
        {
            get
            {
                return _this ?? (_this = (CSARToolsModule)FrameworkApplication.FindModule("Combat_Search_and_Rescue_Tools_Module"));
            }
        }

        #region Overrides
        /// <summary>
        /// Called by Framework when ArcGIS Pro is closing
        /// </summary>
        /// <returns>False to prevent Pro from closing, otherwise True</returns>
        protected override bool CanUnload()
        {
            //TODO - add your business logic
            //return false to ~cancel~ Application close
            return true;
        }

        protected override bool Initialize()
        {
            if (!base.Initialize())
                return false;

            CopyLayoutTemplate();

            return true;
        }
        #endregion Overrides

        private void CopyLayoutTemplate()
        {
            try
            {
                string _templatePath = ApplicationOptions.LayoutOptions.LayoutTemplatePath;
                string _layoutPath = Path.Combine(_templatePath, "Rescue GRG.pagx");
                if (!File.Exists(_templatePath))
                {
                    File.Copy(Path.Combine(_installPath, @"layout_templates\Rescue GRG.pagx"), _layoutPath);
                }
            }
            catch { }
        }
    }
}
