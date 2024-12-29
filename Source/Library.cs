using System.Runtime.InteropServices;

namespace LeidenCD;

// TODO: export the algorithm as a shared library
public static partial class Library
{
    [UnmanagedCallersOnly(EntryPoint = "NETCD_Init")]
    public static void Init()
    {
        // TODO: initialize the algorithm
    }

    // TODO: architecture-agnostic methods
}